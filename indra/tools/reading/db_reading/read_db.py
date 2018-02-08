"""This module provides essential tools to run reading using indra's own
database. This may also be run as a script; for details run:
`python read_pmids_db --help`
"""
from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import logging
import random
import zlib
import pickle
import json
from math import log10, floor, ceil
from datetime import datetime
from indra.tools.reading.util.script_tools import get_parser, make_statements, \
                                             StatementData

logger = logging.getLogger('make_db_readings')
if __name__ == '__main__':
    parser = get_parser(
        'A tool to read and process content from the database.',
        ('A file containing a list of ids of the form <id_type>:<id>. '
         'Note that besided the obvious id types (pmid, pmcid, doi, etc.), '
         'you may use trid and tcid to indicate text ref and text content '
         'ids, respectively. Note that these are specific to the database, '
         'and should thus be used with care.')
        )
    parser.add_argument(
        '-m', '--mode',
        choices=['all', 'unread', 'none'],
        default='unread',
        help=('Set the reading mode. If \'all\', read everything, if '
              '\'unread\', only read content that does not have pre-existing '
              'readings of the same reader and version, if \'none\', only '
              'use pre-existing readings. Default is \'unread\'.')
        )
    parser.add_argument(
        '-t', '--temp',
        default='.',
        help='Select the location of the temp file.'
        )
    parser.add_argument(
        '-o', '--output',
        dest='name',
        help=('Pickle all results and save in files labelled as '
              '<NAME>_<output_type>.pkl.'),
        default=None
        )
    parser.add_argument(
        '-b', '--inner_batch',
        dest='b_in',
        help=('Choose the size of the inner batches, which is the number of '
              'text content entires loaded at a given time, and the number of '
              'entries that are read at a time by a reader. The default is '
              '1,000.'),
        default=1000,
        type=int
        )
    parser.add_argument(
        '-B', '--outer_batch',
        dest='b_out',
        default=10000,
        type=int,
        help=('Select the number of ids to read per outer level batch. This '
              'determines the number of readings/statements uploaded/pickled '
              'at a time, and thus also limits the amount of RAM that will be '
              'used. A larger outer batch means more RAM. The default is '
              '10,000.')
        )
    parser.add_argument(
        '--no_reading_upload',
        help='Choose not to upload the reading output to the database.',
        action='store_true'
        )
    parser.add_argument(
        '--no_statement_upload',
        help='Choose not to upload the statements to the databse.',
        action='store_true'
        )
    parser.add_argument(
        '--force_fulltext',
        help='Make the reader only read full text from the database.',
        action='store_true'
        )
    args = parser.parse_args()
    if args.debug and not args.quiet:
        logger.setLevel(logging.DEBUG)

from indra.db import get_primary_db, formats, texttypes
from indra.db import sql_expressions as sql
from indra.db.util import insert_agents
from indra.tools.reading.readers import get_readers, ReadingData, _get_dir


class ReadDBError(Exception):
    pass


# =============================================================================
# Useful functions
# =============================================================================


def _convert_id_entry(id_entry, allowed_types=None):
    ret = [s.strip() for s in id_entry.split(':')]
    if len(ret) != 2:
        raise ReadDBError('Improperly formatted id entry: \"%s\"' % id_entry)
    ret[0] = ret[0].lower()
    if allowed_types is not None and ret[0] not in allowed_types:
        raise ReadDBError('Invalid id type: \"%s\"' % ret[0])
    return ret


def _enrich_reading_data(reading_data_iter, db=None):
    """Get db ids for all ReadingData objects that correspond to a db ref.

    Note that the objects are modified IN PLACE, so nothing is returned, and if
    a copy of the objects is passed as an argument, this function will have no
    effect. This does nothing if the readings are not in the database.
    """
    logger.debug("Enriching the reading data with database refs.")
    if db is None:
        db = get_primary_db()
    possible_matches = db.select_all(
        'readings',
        db.Readings.text_content_id.in_([rd.tcid for rd in reading_data_iter
                                         if rd.reading_id is None])
        )
    for rdata in reading_data_iter:
        for reading in possible_matches:
            if rdata.matches(reading):
                rdata.reading_id = reading.id
                break
    return


# =============================================================================
# Content Retrieval
# =============================================================================


def get_id_dict(id_str_list):
    """Parse the list of id string into a dict."""
    id_types = get_primary_db().TextRef.__table__.columns.keys()
    id_types.remove('id')
    id_types += ['trid', 'tcid']
    id_dict = {id_type: [] for id_type in id_types}
    for id_entry in id_str_list:
        id_type, id_val = _convert_id_entry(id_entry, id_types)
        if id_type in ['trid', 'tcid']:
            id_dict[id_type].append(int(id_val))
        else:
            id_dict[id_type].append(id_val)
    return id_dict


def get_clauses(id_dict, db):
    """Get a list of clauses to be passed to a db query.

    Note that an empty condition will be returned if id_dict has no ids in it
    (either the dict is empty or all the lists within the dict are empty),
    which will in general have the unexpected effect of selecting everything,
    rather than nothing.

    Parameters
    ----------
    id_dict : dict {id_type: [int or str]}
        A dictionary indexed by the type of id, containing lists of id's of
        that the respective type. If all the lists are empty, or the dict is
        empty, returns an empty condition. Note that id types of 'trid' and
        'tcid' will be mapped to text ref ids and text content ids,
        respectively.
    db : indra.db.DatabaseManager instance
        This instance is only used for forming the query, and will not be
        accessed or queried.

    Returns
    -------
    clause_list : list [sqlalchemy clauses]
        A list of sqlalchemy clauses to be used in query in the form:
        `db.filter_query(<table>, <other clauses>, *clause_list)`.
        If the id_dict has no ids, an effectively empty condition is returned.
    """
    id_condition_list = [getattr(db.TextRef, id_type).in_(id_list)
                         for id_type, id_list in id_dict.items()
                         if len(id_list) and id_type not in ['tcid', 'trid']]
    for id_type, table in [('trid', db.TextRef), ('tcid', db.TextContent)]:
        if id_type in id_dict.keys() and len(id_dict[id_type]):
            int_id_list = [int(i) for i in id_dict[id_type]]
            id_condition_list.append(table.id.in_(int_id_list))
    return [sql.or_(*id_condition_list)]


def get_text_content_summary_string(q, db, num_ids=None):
    """Create a table with some summary data for a query."""
    N_tot = q.count()
    if num_ids is not None:
        logger.info("Found %d text content entires out of %d ids."
                    % (N_tot, num_ids))
    if N_tot > 0:
        log_n = floor(log10(N_tot))
    else:
        log_n = 1
    cols = list(formats.values()) + ['tot']
    col_fmt = ' %%%ds' % max(4, log_n)
    cols_strs = [col_fmt % fmt for fmt in cols]
    ret_str = 'Summary Statisitics:\n' + ' '*10 + ''.join(cols_strs) + '\n'
    col_counts = dict.fromkeys(formats.values())
    col_counts['tot'] = []
    for texttype in texttypes.values():
        line = '%8s: ' % texttype
        counts = []
        for text_format in cols[:-1]:
            if col_counts[text_format] is None:
                col_counts[text_format] = []
            c = q.filter(
                db.TextContent.text_type == texttype,
                db.TextContent.format == text_format
                ).count()
            line += col_fmt % c
            counts.append(c)
            col_counts[text_format].append(c)
        line += col_fmt % sum(counts)
        ret_str += line + '\n'
        col_counts['tot'].append(sum(counts))
    ret_str += '%8s: ' % 'total' + ''.join([col_fmt % sum(col_counts[col])
                                            for col in cols])
    return ret_str


def get_content_query(ids, readers, db=None, force_fulltext=False,
                      force_read=False, debug=False, print_summary=False):
    """Construct a query to access all the content that will be read.

    If ids is not 'all', and does not contain any ids, None is returned.

    Parameters
    ----------
    ids : 'all' or dict {<id type> : [str/int]}
        If 'all', then all the content will be included in the query. Otherwise
        a the content will be constrained to that corresponding to the ids in
        id_dict, which are matched using text refs.
    readers : list [Reader child instances]
        A list of the reader objects, which contain the required metadata (name
        and version of the reader) used to find content that needs to be read.
    db : indra.db.DatabaseManager instance
        Optional, default None, in which case the primary database is used. If
        specified, the alternative database will be used. This function should
        not alter the database.
    force_fulltext : bool
        Optional, default False - If True, only fulltext content will be read,
        as opposed to including abstracts.
    force_read : bool
        Optional, default False - If True, all content will be returned,
        whether it has been read or not.

    Returns
    -------
    tc_tbr_query : sqlalchemy query object or None
        The query of the text content to be read (tc_tbr). If there are no ids
        contained in ids, or it is not 'all', return None.
    """
    if debug:
        logger.setLevel(logging.DEBUG)
    if db is None:
        db = get_primary_db()
    logger.debug("Got db handle.")

    # These allow conditions on different tables to equal conditions on the
    # dependent tables.
    tc_tr_binding = db.TextContent.text_ref_id == db.TextRef.id
    rd_tc_binding = db.Readings.text_content_id == db.TextContent.id

    # Begin the list of clauses with the binding between text content and
    # text refs.
    clauses = [tc_tr_binding]

    # Add a fulltext requirement, if applicable.
    if force_fulltext:
        clauses.append(db.TextContent.text_type == texttypes.FULLTEXT)

    # If we are actually getting anything, else we return None.
    if ids == 'all' or any([len(id_list) > 0 for id_list in ids.values()]):
        if ids is not 'all':
            clauses += get_clauses(ids, db)

        # Get the text content query object
        tc_query = db.filter_query(
            db.TextContent,
            *clauses
            ).distinct()

        if not force_read:
            logger.debug("Getting content to be read.")
            # Each sub query is a set of content that has been read by one of
            # the readers.
            tc_q_subs = [tc_query.filter(rd_tc_binding, r.matches_clause(db))
                         for r in readers]
            tc_tbr_query = tc_query.except_(sql.intersect(*tc_q_subs))
        else:
            logger.debug('All content will be read (force_read).')
            tc_tbr_query = tc_query

        if print_summary:
            try:
                logger.debug("Going to try to make a nice summary...")
                logger.info(get_text_content_summary_string(tc_tbr_query, db))
            except Exception:
                logger.debug("Could not print summary of results.")
    else:
        logger.debug("No ids in id_dict, so no query formed.")
        return None

    return tc_tbr_query.distinct()


def get_readings_query(ids, readers, db=None, force_fulltext=False):
    """Create a query to access all the relevant existing readings.

    Note that if ids is not 'all' and ids is a dict with no ids in it,
    this function returns None.

    Parameters
    ----------
    ids : 'all' or dict {<id_type> : [str/int]}
        If 'all', then all possible readings in the database matching the given
        readers and other conditions will be returned. Otherwise, only those
        that correspond to one of the ids in ids dict will be contained. If an
        ids dict has no ids in it, None is returned.
    readers : list [Reader child instances]
        A list of the readers whose names and versions you wish to match in the
        readings queried from the database.
    db : indra.db.DatabaseManager instance
        Optional, default None, in which case the primary database is used. If
        specified, the alternative database will be used. This function should
        not alter the database.
    force_fulltext : bool
        Optional, default False - If True, only readings corresponding to
        fulltext content will be read, as opposed to including readings created
        from abstracts.

    Returns
    -------
    readings_query : sql query instance or None
        Returns a query that can be used to access the specified content, or
        else None if no content was specified.
    """
    if db is None:
        db = get_primary_db()
    clauses = [
        # Bind conditions on readings to conditions on content.
        db.Readings.text_content_id == db.TextContent.id,

        # Bind text content to text refs
        db.TextContent.text_ref_id == db.TextRef.id,

        # Check if at least one of the readers has read the content
        sql.or_(*[reader.matches_clause(db) for reader in readers])
        ]
    if force_fulltext:
        clauses.append(db.TextContent.text_type == texttypes.FULLTEXT)

    if ids == 'all' or any([id_list for id_list in ids.values()]):
        if ids != 'all':
            clauses += get_clauses(ids, db)

        readings_query = db.filter_query(
            db.Readings,

            # Bind conditions on readings to conditions on content.
            db.Readings.text_content_id == db.TextContent.id,

            # Bind text content to text refs
            db.TextContent.text_ref_id == db.TextRef.id,

            # Check if at least one of the readers has read the content
            sql.or_(*[reader.matches_clause(db) for reader in readers]),

            # Conditions generated from the list of ids. These include a
            # text-ref text-content binding to connect with id data.
            *clauses
            )
    else:
        return None

    return readings_query.distinct()


# =============================================================================
# Core Reading Functions
# =============================================================================


def make_db_readings(id_dict, readers, batch_size=1000, force_fulltext=False,
                     force_read=False, skip_dict=None, db=None, **kwargs):
    """Read contents retrieved from the database.

    The content will be retrieved in batchs, given by the `batch` argument.
    This prevents the system RAM from being overloaded.

    Parameters
    ----------
    id_dict : dict {<id_type>:[<id value>, ...]}
        A dict of lists of the id's to be read, keyed by id_type.
    readers : list of reader objects
        A list of the readers that will be use, for example ['reach'] if you
        wanted to use the reach reader.
    batch_size : int
        The number of content entries read for each batch. Default 1000.
    force_fulltext : bool
        If True, only get fulltext content from the database. Default False.
    force_read : bool
        If True, read even if text_content id is found in skip_dict.
    skip_dict : dict {<reader> : list [int]}
        A dict containing text content id's to be skipped.
    db : indra.db.DatabaseManager instance
        A handle to a database. Default None; if None, a handle to the primary
        database (see indra.db) is retrieved.

    Other keyword arguments are passed to the `read` methods of the readers.

    Returns
    -------
    outputs : list of ReadingData instances
        The results of the readings with relevant metadata.
    """
    if db is None:
        db = get_primary_db()

    # Get the iterator.
    logger.debug("Getting iterator.")
    tc_read_q = get_content_query(
        id_dict,
        readers,
        db=db,
        force_fulltext=force_fulltext,
        force_read=force_read
        )
    logger.debug("Begginning to iterate.")
    batch_list_dict = {r.name: [] for r in readers}
    new_outputs = []
    if tc_read_q is not None:
        for text_content in tc_read_q.yield_per(batch_size):
            # The get_content function returns an iterator which yields
            # results in batches, so as not to overwhelm RAM. We need to read
            # in batches for much the same reaason.
            for r in readers:
                if not force_read:
                    if skip_dict is not None:
                        if text_content.id in skip_dict[r.name]:
                            continue
                    else:
                        # Try to get a previous reading from this reader.
                        reading = db.select_one(
                            db.Readings,
                            db.Readings.text_content_id == text_content.id,
                            r.matches_clause(db)
                            )
                        if reading is not None:
                            continue
                batch_list_dict[r.name].append(text_content)

                if (len(batch_list_dict[r.name])+1) % batch_size is 0:
                    # TODO: this is a bit cludgy...maybe do this better?
                    # Perhaps refactor read_content.
                    logger.debug("Reading batch of files for %s." % r.name)
                    results = r.read(batch_list_dict[r.name], **kwargs)
                    if results is not None:
                        new_outputs += results
                    batch_list_dict[r.name] = []
        logger.debug("Finished iteration.")
        # Pick up any stragglers.
        for r in readers:
            if len(batch_list_dict[r.name]) > 0:
                logger.debug("Reading remaining files for %s." % r.name)
                results = r.read(batch_list_dict[r.name], **kwargs)
                if results is not None:
                    new_outputs += results
    return new_outputs


def get_db_readings(id_dict, readers, force_fulltext=False, batch_size=1000,
                    db=None):
    """Get readings from the database."""
    if db is None:
        db = get_primary_db()

    # Get any previous readings. Note that we do this BEFORE posting the new
    # readings. Otherwise we would have duplicates.
    previous_readings_query = get_readings_query(
        id_dict,
        readers,
        db=db,
        force_fulltext=force_fulltext
        )
    if previous_readings_query is not None:
        prev_readings = [
            ReadingData(
                r.text_content_id,
                r.reader,
                r.reader_version,
                r.format,
                json.loads(zlib.decompress(r.bytes, 16+zlib.MAX_WBITS)
                           .decode('utf8')),
                r.id
                )
            for r in previous_readings_query.yield_per(batch_size)
            ]
    else:
        prev_readings = []
    return prev_readings


def upload_readings(output_list, db=None):
    """Put the reading output on the database."""
    if db is None:
        db = get_primary_db()

    # Create the list of records to be copied, ensuring no uniqueness conflicts
    r_list = db.select_all(
        db.Readings,
        db.Readings.text_content_id.in_([rd.tcid for rd in output_list])
        )
    exisiting_tcid_set = set([r.text_content_id for r in r_list])
    upload_list = []
    for reading_data in output_list:
        # First check if this tcid is even in the set of existing tcids in the
        # readings table.
        if reading_data.tcid in exisiting_tcid_set:
            r_tcid_list = [r for r in r_list
                           if r.text_content_id == reading_data.tcid]
            # Now check for any exact matches:
            if any([reading_data.matches(r) for r in r_tcid_list]):
                continue

        # If there were no conflicts, we can add this to the copy list.
        upload_list.append(reading_data.make_tuple())

    # Copy into the database.
    logger.info("Adding %d/%d reading entries to the database." %
                (len(upload_list), len(output_list)))
    db.copy('readings', upload_list, ReadingData.get_cols())
    return


def produce_readings(id_dict, reader_list, verbose=False, read_mode='unread',
                     force_fulltext=False, batch_size=1000, no_upload=False,
                     pickle_file=None, db=None, log_readers=True):
    """Produce the reading output for the given ids, and upload them to db.

    This function will also retrieve pre-existing readings from the database,
    thus improving performance.

    Parameters
    ----------
    id_dict : dict {<id_type>:[<id value>, ...]}
        A dict of lists of the id's to be read, keyed by id_type.
    reader_list : list [Reader]
        A list of Reader descendents to be used in reading.
    verbose : bool
        Optional, default False - If True, log and print the output of the
        commandline reader utilities, if False, don't.
    read_mode : str : 'all', 'unread', or 'none'
        Optional, default 'undread' - If 'all', read everything (generally
        slow); if 'unread', only read things that were undread (as fast as you
        can be while still getting everything); if 'none', don't read, and only
        get existing readings.
    force_fulltext : bool
        Optional, default False - If True, only read fulltext article, ignoring
        abstracts.
    batch_size : int
        Optional, default 1000 - The number of text content entries to be
        yielded by the database at a given time.
    no_read : bool
        Optional, default False - If True, do not perform any new readings, and
        only retrieve existing readings from the database.
    no_upload : bool
        Optional, default False - If True, do not upload content to the
        database.
    pickle_file : str or None
        Optional, default None - otherwise the path to a file in which the
        reading data will be saved.
    db : indra.db.DatabaseManager instance
        Optional, default is None, in which case the primary database provided
        by `get_primary_db` function is used. Used to interface with a
        different databse.

    Returns
    -------
    outputs : list [ReadingData]
        A list of the outputs of the readings in the form of ReadingData
        instances.
    """
    logger.debug("Producing readings in %s mode." % read_mode)
    if db is None:
        db = get_primary_db()

    prev_readings = []
    skip_reader_tcid_dict = None
    if read_mode != 'all':
        prev_readings = get_db_readings(id_dict, reader_list, force_fulltext,
                                        batch_size, db=db)
        skip_reader_tcid_dict = {r.name: [] for r in reader_list}
        logger.info("Found %d pre-existing readings." % len(prev_readings))
        if read_mode != 'none':
            for rd in prev_readings:
                skip_reader_tcid_dict[rd.reader].append(rd.tcid)
    outputs = []
    if read_mode != 'none':
        outputs = make_db_readings(id_dict, reader_list, verbose=verbose,
                                   skip_dict=skip_reader_tcid_dict, db=db,
                                   force_fulltext=force_fulltext,
                                   force_read=(read_mode == 'all'),
                                   batch_size=batch_size, log=log_readers)
        logger.info("Made %d new readings." % len(outputs))

    if not no_upload:
        try:
            upload_readings(outputs, db=db)
        except Exception as e:
            logger.exception(e)
            if pickle_file is None:
                pickle_file = ("failure_reading_dump_%s.pkl" 
                               % datetime.now().strftime('%Y%m%d_%H%M%S'))
            logger.error("Cound not upload readings. Results are pickled in: "
                         + pickle_file)

    outputs += prev_readings

    if pickle_file is not None:
        with open(pickle_file, 'wb') as f:
            pickle.dump([output.make_tuple() for output in outputs], f)
        print("Reading outputs stored in %s." % pickle_file)

    return outputs


# =============================================================================
# Statement Processing
# =============================================================================


def upload_statements(stmt_data_list, db=None):
    """Upload the statements to the database."""
    if db is None:
        db = get_primary_db()

    logger.info("Uploading %d statements to the database." %
                len(stmt_data_list))
    db.copy('statements', [s.make_tuple() for s in stmt_data_list],
            StatementData.get_cols())

    logger.info("Uploading agents to the database.")
    reading_id_set = set([sd.reading_id for sd in stmt_data_list])
    if len(reading_id_set):
        insert_agents(db, [sd.statement for sd in stmt_data_list],
                         db.Statements.reader_ref.in_(reading_id_set))
    return


def produce_statements(output_list, enrich=True, no_upload=False,
                       pickle_file=None, n_proc=1, db=None):
    """Convert the reader output into a list of StatementData instances."""
    if db is None:
        db = get_primary_db()

    if enrich:
        _enrich_reading_data(output_list, db=db)

    stmt_data_list = make_statements(output_list, n_proc)

    if not no_upload:
        try:
            upload_statements(stmt_data_list, db=db)
        except Exception as e:
            logger.exception(e)
            if pickle_file is None:
                pickle_file = ("failure_stmt_dump_%s.pkl"
                               % datetime.now().strftime('%Y%m%d_%H%M%S'))
            logger.error("Could not upload statements. Results pickled in: %s."
                         % pickle_file)
    if pickle_file is not None:
        with open(pickle_file, 'wb') as f:
            pickle.dump([sd.statement for sd in stmt_data_list], f)
        print("Statements pickled in %s." % pickle_file)

    return stmt_data_list


# =============================================================================
# Main for script use
# =============================================================================


if __name__ == "__main__":
    # Process the arguments. =================================================

    # Get the ids.
    with open(args.input_file, 'r') as f:
        input_lines = f.readlines()
    logger.info("Found %d ids." % len(input_lines))

    # Select only a sample of the lines, if sample is chosen.
    if args.n_samp is not None:
        input_lines = random.sample(input_lines, args.n_samp)
    else:
        random.shuffle(input_lines)

    # If a range is specified, only use that range.
    if args.range_str is not None:
        start_idx, end_idx = [int(n) for n in args.range_str.split(':')]
        input_lines = input_lines[start_idx:end_idx]

    # Get the outer batch.
    B = args.b_out
    n_max = int(ceil(float(len(input_lines))/B))

    # Create a single base directory
    base_dir = _get_dir(args.temp, 'run_%s' % ('_and_'.join(args.readers)))

    # Get the readers objects.
    readers = [reader_class(base_dir=base_dir, n_proc=args.n_proc)
               for reader_class in get_readers()
               if reader_class.name.lower() in args.readers]

    # Set the verbosity. The quiet argument overrides the verbose argument.
    verbose = args.verbose and not args.quiet

    for n in range(n_max):
        logger.info("Beginning outer batch %d/%d. ------------" % (n+1, n_max))

        # Get the pickle file names.
        if args.name is not None:
            reading_pickle = args.name + '_readings_%d.pkl' % n
            stmts_pickle = args.name + '_stmts_%d.pkl' % n
        else:
            reading_pickle = None
            stmts_pickle = None

        # Get the dict of ids.
        id_dict = get_id_dict(input_lines[B*n:B*(n+1)])

        # Read everything ====================================================
        outputs = produce_readings(id_dict, readers, verbose=verbose,
                                   read_mode=args.mode, batch_size=args.b_in,
                                   force_fulltext=args.force_fulltext,
                                   no_upload=args.no_reading_upload,
                                   pickle_file=reading_pickle)

        # Convert the outputs to statements ==================================
        produce_statements(outputs, no_upload=args.no_statement_upload,
                           pickle_file=stmts_pickle, n_proc=args.n_proc)
