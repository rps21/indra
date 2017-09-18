from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import boto3
import logging
import botocore.session
from time import sleep
from indra.literature import elsevier_client as ec

bucket_name = 'bigmech'

logger = logging.getLogger('aws_reading')

def wait_for_complete(job_list, poll_interval=10):
    """Return when all jobs in the given list finished.

    If not job list is given, return when all jobs in queue finished.

    Parameters
    ----------
    job_list : Optional[list(dict)]
        A list of jobID-s in a dict, as returned by the submit function.
        Example: [{'jobId': 'e6b00f24-a466-4a72-b735-d205e29117b4'}, ...]
        If not given, this function will return if all jobs completed.
    poll_interval : Optional[int]
        The time delay between API calls to check the job statuses.
    """
    def get_jobs_by_status(status, job_filter=None):
        res = batch_client.list_jobs(jobQueue='run_reach_queue',
                                     jobStatus=status)
        ids = [job['jobId'] for job in res['jobSummaryList']]
        if job_filter:
            ids = [job_id for job_id in ids if job_id in job_filter]
        return ids

    job_list = [job['jobId'] for job in job_list]

    batch_client = boto3.client('batch')

    total_time = 0
    while True:
        not_done = []
        for status in ('SUBMITTED', 'PENDING', 'RUNNABLE', 'STARTING',
                       'RUNNING'):
            not_done += get_jobs_by_status(status, job_list)
        failed = get_jobs_by_status('FAILED', job_list)
        done = get_jobs_by_status('SUCCEEDED', job_list)

        logger.info('(%d s)=(not done: %d, failed: %d, done: %d)' %
              (total_time, len(not_done), len(failed), len(done)))

        if job_list:
            if (len(failed) + len(done)) == len(job_list):
                return 0
        else:
            if (len(failed) + len(done) > 0) and (len(not_done) == 0):
                return 0

        sleep(poll_interval)
        total_time += poll_interval

def get_environment():
    # Get AWS credentials
    # http://stackoverflow.com/questions/36287720/boto3-get-credentials-dynamically
    session = botocore.session.get_session()
    access_key = session.get_credentials().access_key
    secret_key = session.get_credentials().secret_key

    # Get the Elsevier keys from the Elsevier client
    environment_vars = [
            {'name': ec.api_key_env_name,
             'value': ec.elsevier_keys.get('X-ELS-APIKey')},
            {'name': ec.inst_key_env_name,
             'value': ec.elsevier_keys.get('X-ELS-Insttoken')},
            {'name': 'AWS_ACCESS_KEY_ID',
             'value': access_key},
            {'name': 'AWS_SECRET_ACCESS_KEY',
             'value': secret_key}]
    return environment_vars


def submit_run_reach(basename, pmid_list_filename, start_ix=None, end_ix=None,
                     pmids_per_job=3000, num_tries=2):
    # Upload the pmid_list to Amazon S3
    pmid_list_key = 'reading_results/%s/pmids' % basename
    s3_client = boto3.client('s3')
    s3_client.upload_file(pmid_list_filename, 'bigmech', pmid_list_key)

    # If no end index is specified, read all the PMIDs
    if end_ix is None:
        with open(pmid_list_filename, 'rt') as f:
            lines = f.readlines()
            end_ix = len(lines)

    if start_ix is None:
        start_ix = 0

    # Get environment variables
    environment_vars = get_environment()

    # Iterate over the list of PMIDs and submit the job in chunks
    batch_client = boto3.client('batch')
    job_list = []
    for job_start_ix in range(start_ix, end_ix, pmids_per_job):
        job_end_ix = job_start_ix + pmids_per_job
        if job_end_ix > end_ix:
            job_end_ix = end_ix
        job_name = '%s_%d_%d' % (basename, job_start_ix, job_end_ix)
        command_list = ['python', '-m',
                        'indra.tools.reading.run_reach_on_pmids_aws',
                        basename, '/tmp', '16', str(job_start_ix),
                        str(job_end_ix)]
        print(command_list)
        job_info = batch_client.submit_job(jobName=job_name,
                       jobQueue='run_reach_queue',
                       jobDefinition='run_reach_jobdef',
                       containerOverrides={'environment': environment_vars,
                                           'command': command_list})#,
                       #retryStrategy={'attempts': num_tries})
        job_list.append({'jobId': job_info['jobId']})
    return job_list


def submit_combine(basename, job_ids=None):
    # Get environment variables
    environment_vars = get_environment()

    job_name = '%s_combine_reach' % basename
    command_list = ['python', '-m',
                    'indra.tools.reading.assemble_reach_stmts_aws',
                    basename]
    print(command_list)
    kwargs = {'jobName': job_name, 'jobQueue': 'run_reach_queue',
              'jobDefinition': 'run_reach_jobdef',
              'containerOverrides': {'environment': environment_vars,
                                     'command': command_list,
                                     'memory': 60000, 'vcpus': 1}}
    if job_ids:
        kwargs['dependsOn'] = job_ids
    batch_client = boto3.client('batch')
    batch_client.submit_job(**kwargs)


if __name__ == '__main__':
    import sys
    import argparse

    # Create the top-level parser
    parser = argparse.ArgumentParser('submit_reading_pipeline_aws.py',
            description='Run machine reading with REACH using AWS Batch.')
    subparsers = parser.add_subparsers(title='Job Type',
            help='Type of jobs to submit.')
    subparsers.required = True
    subparsers.dest = 'job_type'
    parent_read_parser = argparse.ArgumentParser(add_help=False)
    parent_read_parser.add_argument('basename',
        help='Defines job names and S3 keys')
    parent_read_parser.add_argument('pmid_file',
        help='Path to file containing PMIDs to read')
    parent_read_parser.add_argument('--start_ix', type=int,
        help='Start index of PMIDs to read.')
    parent_read_parser.add_argument('--end_ix', type=int,
        help='End index of PMIDs to read (default: read all PMIDs)')
    parent_read_parser.add_argument('--force_read', action='store_true',
        help='Read papers even if previously read by current REACH.')
    parent_read_parser.add_argument('--force_fulltext', action='store_true',
        help='Get full text content even if content already on S3.')
    parent_read_parser.add_argument('--pmids_per_job', default=3000, type=int,
        help='Number of PMIDs to read for each AWS Batch job.')
    parent_read_parser.add_argument('--num_tries', default=2, type=int,
        help='Maximum number of times to try running job.')
    read_parser = subparsers.add_parser('read', parents=[parent_read_parser],
        help='Run REACH and cache INDRA Statements on S3.',
        description='Run REACH and cache INDRA Statements on S3.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    combine_parser = subparsers.add_parser('combine',
        help='Combine INDRA Statement subsets into a single file.',
        description='Combine INDRA Statement subsets into a single file.')
    combine_parser.add_argument('basename',
        help='Defines job name and S3 keys')
    full_parser = subparsers.add_parser('full',
        parents=[parent_read_parser],
        help='Run REACH and combine INDRA Statements when done.',
        description='Run REACH and combine INDRA Statements when done.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    args = parser.parse_args()

    if args.job_type == 'read':
        job_ids = submit_run_reach(args.basename, args.pmid_file,
                                   args.start_ix, args.end_ix,
                                   args.pmids_per_job)
    elif args.job_type == 'combine':
        submit_combine(args.basename)
    elif args.job_type == 'full':
        job_ids = submit_run_reach(args.basename, args.pmid_file,
                                   args.start_ix, args.end_ix,
                                   args.pmids_per_job)
        submit_combine(args.basename, job_ids)
    else:
        print('job_type must be one of ("read", "combine", "full")')
        sys.exit(1)


