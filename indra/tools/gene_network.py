from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import pickle
import logging
from indra.sources import bel, biopax
import indra.tools.assemble_corpus as ac
from indra.preassembler import Preassembler
from indra.preassembler.hierarchy_manager import hierarchies
from indra.preassembler.sitemapper import default_mapper as sm

logger = logging.getLogger('gene_network')

class GeneNetwork(object):
    """Build a set of INDRA statements for a given gene list from databases.

    Parameters
    ----------
    gene_list : string
        List of gene names.
    basename : string or None (default)
        Filename prefix to be used for caching of intermediates (Biopax OWL
        file, pickled statement lists, etc.). If None, no results are cached
        and no cached files are used.

    Attributes
    ----------
    gene_list : string
        List of gene names
    basename : string or None
        Filename prefix for cached intermediates, or None if no cached used.
    results : dict
        Dict containing results of preassembly (see return type for
        :py:meth:`run_preassembly`.
    """
    def __init__(self, gene_list, basename=None):
        if not gene_list:
            raise ValueError("Gene list must contain at least one element.")
        self.gene_list = gene_list
        self.basename = basename

    def get_bel_stmts(self, filter=False):
        """Get relevant statements from the BEL large corpus.

        Performs a series of neighborhood queries and then takes the union of
        all the statements. Because the query process can take a long time for
        large gene lists, the resulting list of statements are cached in a
        pickle file with the filename `<basename>_bel_stmts.pkl`.  If the
        pickle file is present, it is used by default; if not present, the
        queries are performed and the results are cached.

        Parameters
        ----------
        filter : bool
            If True, includes only those statements that exclusively mention
            genes in :py:attr:`gene_list`. Default is False. Note that the
            full (unfiltered) set of statements are cached.

        Returns
        -------
        list of :py:class:`indra.statements.Statement`
            List of INDRA statements extracted from the BEL large corpus.
        """
        if self.basename is not None:
            bel_stmt_path = '%s_bel_stmts.pkl' % self.basename
        # Check for cached BEL stmt file
        if self.basename is not None and os.path.isfile(bel_stmt_path):
            logger.info("Loading BEL statements from %s" % bel_stmt_path)
            with open(bel_stmt_path, 'rb') as f:
                bel_statements = pickle.load(f)
        # No cache, so perform the queries
        else:
            bel_statements = []
            for gene in self.gene_list:
                logger.info("Getting BEL statements for gene %s" % gene)
                bel_proc = bel.process_ndex_neighborhood([gene])
                if bel_proc is not None:
                    bel_statements += bel_proc.statements
            # Save to pickle file if we're caching
            if self.basename is not None:
                with open(bel_stmt_path, 'wb') as f:
                    pickle.dump(bel_statements, f, protocol=2)
        # Optionally filter out statements not involving only our gene set
        if filter:
            if len(self.gene_list) > 1:
                bel_statements = ac.filter_gene_list(bel_statements,
                                                     self.gene_list, 'one')
            else:
                bel_statements = ac.filter_gene_list(bel_statements,
                                                     self.gene_list, 'all')
        return bel_statements

    def get_biopax_stmts(self, filter=False, query='pathsbetween',
                         database_filter=None):
        """Get relevant statements from Pathway Commons.

        Performs a "paths between" query for the genes in :py:attr:`gene_list`
        and uses the results to build statements. This function caches two
        files: the list of statements built from the query, which is cached in
        `<basename>_biopax_stmts.pkl`, and the OWL file returned by the Pathway
        Commons Web API, which is cached in `<basename>_pc_pathsbetween.owl`.
        If these cached files are found, then the results are returned based
        on the cached file and Pathway Commons is not queried again.

        Parameters
        ----------
        filter : Optional[bool]
            If True, includes only those statements that exclusively mention
            genes in :py:attr:`gene_list`. Default is False.
        query : Optional[str]
            Defined what type of query is executed. The two options are
            'pathsbetween' which finds paths between the given list of genes
            and only works if more than 1 gene is given, and 'neighborhood'
            which searches the immediate neighborhood of each given gene.
            Note that for pathsbetween queries with more thatn 60 genes, the
            query will be executed in multiple blocks for scalability.
        database_filter: Optional[list[str]]
            A list of PathwayCommons databases to include in the query.

        Returns
        -------
        list of :py:class:`indra.statements.Statement`
            List of INDRA statements extracted from Pathway Commons.
        """
        # If we're using a cache, initialize the appropriate filenames
        if self.basename is not None:
            biopax_stmt_path = '%s_biopax_stmts.pkl' % self.basename
            biopax_ras_owl_path = '%s_pc_pathsbetween.owl' % self.basename
        # Check for cached Biopax stmt file at the given path
        # if it's there, return the statements from the cache
        if self.basename is not None and os.path.isfile(biopax_stmt_path):
            logger.info("Loading Biopax statements from %s" % biopax_stmt_path)
            with open(biopax_stmt_path, 'rb') as f:
                bp_statements = pickle.load(f)
            return bp_statements
        # Check for cached file before querying Pathway Commons Web API
        if self.basename is not None and os.path.isfile(biopax_ras_owl_path):
            logger.info("Loading Biopax from OWL file %s" % biopax_ras_owl_path)
            bp = biopax.process_owl(biopax_ras_owl_path)
        # OWL file not found; do query and save to file
        else:
            if (len(self.gene_list) < 2) and (query == 'pathsbetween'):
                logger.warning('Using neighborhood query for one gene.')
                query = 'neighborhood'
            if query == 'pathsbetween':
                if len(self.gene_list) > 60:
                    block_size = 60
                else:
                    block_size = None
                bp = biopax.process_pc_pathsbetween(self.gene_list,
                                                database_filter=database_filter,
                                                block_size=block_size)
            elif query == 'neighborhood':
                bp = biopax.process_pc_neighborhood(self.gene_list,
                                                database_filter=database_filter)
            else:
                logger.error('Invalid query type: %s' % query)
                return []
            # Save the file if we're caching
            if self.basename is not None:
                bp.save_model(biopax_ras_owl_path)
        # Save statements to pickle file if we're caching
        if self.basename is not None:
            with open(biopax_stmt_path, 'wb') as f:
                pickle.dump(bp.statements, f, protocol=2)
        # Optionally filter out statements not involving only our gene set
        if filter:
            if len(self.gene_list) > 1:
                bp_statements = ac.filter_gene_list(bp_statements,
                                                     self.gene_list, 'one')
            else:
                bp_statements = ac.filter_gene_list(bp_statements,
                                                     self.gene_list, 'all')
        return bp.statements

    def get_statements(self, filter=False):
        """Return the combined list of statements from BEL and Pathway Commons.

        Internally calls :py:meth:`get_biopax_stmts` and
        :py:meth:`get_bel_stmts`.

        Parameters
        ----------
        filter : bool
            If True, includes only those statements that exclusively mention
            genes in :py:attr:`gene_list`. Default is False.

        Returns
        -------
        list of :py:class:`indra.statements.Statement`
            List of INDRA statements extracted the BEL large corpus and Pathway
            Commons.
        """
        bp_stmts = self.get_biopax_stmts(filter=filter)
        bel_stmts = self.get_bel_stmts(filter=filter)

        return bp_stmts + bel_stmts

    def run_preassembly(self, stmts, print_summary=True):
        """Run complete preassembly procedure on the given statements.

        Results are returned as a dict and stored in the attribute
        :py:attr:`results`. They are also saved in the pickle file
        `<basename>_results.pkl`.

        Parameters
        ----------
        stmts : list of :py:class:`indra.statements.Statement`
            Statements to preassemble.
        print_summary : bool
            If True (default), prints a summary of the preassembly process to
            the console.

        Returns
        -------
        dict
            A dict containing the following entries:

            - `raw`: the starting set of statements before preassembly.
            - `duplicates1`: statements after initial de-duplication.
            - `valid`: statements found to have valid modification sites.
            - `mapped`: mapped statements (list of
              :py:class:`indra.preassembler.sitemapper.MappedStatement`).
            - `mapped_stmts`: combined list of valid statements and statements
              after mapping.
            - `duplicates2`: statements resulting from de-duplication of the
              statements in `mapped_stmts`.
            - `related2`: top-level statements after combining the statements
              in `duplicates2`.
        """
        # First round of preassembly: remove duplicates before sitemapping
        pa1 = Preassembler(hierarchies, stmts)
        logger.info("Combining duplicates")
        pa1.combine_duplicates()
        # Map sites
        logger.info("Mapping sites")
        (valid, mapped) = sm.map_sites(pa1.unique_stmts)
        # Combine valid and successfully mapped statements into single list
        correctly_mapped_stmts = []
        for ms in mapped:
            if all([True if mm[1] is not None else False
                         for mm in ms.mapped_mods]):
                correctly_mapped_stmts.append(ms.mapped_stmt)
        mapped_stmts = valid + correctly_mapped_stmts 
        # Second round of preassembly: de-duplicate and combine related
        pa2 = Preassembler(hierarchies, mapped_stmts)
        logger.info("Combining duplicates again")
        pa2.combine_duplicates()
        pa2.combine_related()
        # Fill out the results dict
        self.results = {}
        self.results['raw'] = stmts
        self.results['duplicates1'] = pa1.unique_stmts
        self.results['valid'] = valid
        self.results['mapped'] = mapped
        self.results['mapped_stmts'] = mapped_stmts
        self.results['duplicates2'] = pa2.unique_stmts
        self.results['related2'] = pa2.related_stmts
        # Print summary
        if print_summary:
            logger.info("\nStarting number of statements: %d" % len(stmts))
            logger.info("After duplicate removal: %d" % len(pa1.unique_stmts))
            logger.info("Unique statements with valid sites: %d" % len(valid))
            logger.info("Unique statements with invalid sites: %d" %
                        len(mapped))
            logger.info("After post-mapping duplicate removal: %d" %
                        len(pa2.unique_stmts))
            logger.info("After combining related statements: %d" %
                        len(pa2.related_stmts))
        # Save the results if we're caching
        if self.basename is not None:
            results_filename = '%s_results.pkl' % self.basename
            with open(results_filename, 'wb') as f:
                pickle.dump(self.results, f, protocol=2)
        return self.results

