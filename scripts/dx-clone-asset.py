#!/usr/bin/env python
from __future__ import print_function
import argparse
import subprocess
import multiprocessing
import time
import functools
import sys
import json

import dxpy

ASSET_PROJECT_PREFIX = 'DNAnexus Assets for'
URL_DURATION = 60 * 60 * 24
SLEEP_TIME = 5
CLONE_ASSET_APP_NAME = '_clone_asset'
CLONE_ASSET_APP = dxpy.find_one_app(zero_ok=False, more_ok=False, name=CLONE_ASSET_APP_NAME, return_handler=True)

# Get the set of supported regions
SUPPORTED_REGIONS = set()
user_description = dxpy.api.user_describe(dxpy.whoami())
if user_description['billTo'].startswith('user-'):
    SUPPORTED_REGIONS = set(user_description['permittedRegions'])
elif user_description['billTo'].startswith('org-'):
    SUPPORTED_REGIONS = set(dxpy.api.org_describe(user_description['billTo'])['permittedRegions'])


def _parse_args():
    """
    Parse the input arguments.
    """
    ap = argparse.ArgumentParser(description='Clone an asset')

    ap.add_argument('--record',
                    help='Record-id of asset to clone.',
                    required=True)
    ap.add_argument('--regions',
                    help='Regions to clone asset into.  Permitted regions are:\n[{0}]'.format(', '.join(SUPPORTED_REGIONS)),
                    choices=SUPPORTED_REGIONS,
                    default=SUPPORTED_REGIONS,
                    nargs='+',
                    metavar='',
                    required=False)
    ap.add_argument('--num-retries',
                    help='Number of attempts to transfer the asset to a given region.',
                    default=0,
                    type=int,
                    required=False)
    ap.add_argument('--priority',
                    help='Priority with which to run the clone_asset app',
                    choices = ['normal', 'high'],
                    required=False)

    return ap.parse_args()


def _find_asset_project(region):
    """
    Returns the asset project for the given region, or None if a problem arrises.
    """
    project_name = '{0} {1}'.format(ASSET_PROJECT_PREFIX, region)

    # Try to find the asset project for the given region.
    # If more than 1 project with the asset project name is found
    # for the given region, or if no project is found and one
    # can't be created, return None indicating there was a problem.
    try:
        cmd = 'dx find projects --level CONTRIBUTE --name "{0}" --region "{1}" --brief '
        projects = subprocess.check_output(cmd.format(project_name, region), shell=True).strip()
        if projects == '':
            cmd = 'dx new project --region "{0}" "{1}" --brief '
            project = subprocess.check_output(cmd.format(region, project_name), shell=True).strip()
        else:
            projects = projects.split('\n')
            return projects[0]
    except:
        pass

    return None


def _clone_asset_into_region(region, record_name, asset_properties, asset_file_name, url, num_retries, q, priority):
    """
    Run the _clone_asset app to clone the given asset information into a new asset in the given region.
    The new asset will live in a project in the given region with a prefix given by ASSET_PROJECT_PREFIX.
    This function will attempt to re-run the transfer app num_retries times before finally giving up.

    The function will return the record_id of the new asset if successful, or None if it is not successful.
    """
    # Get the official asset project for the given region.
    project = _find_asset_project(region)
    # If no official asset project is found and one can't be created,
    # just return None.
    if project is None:
        return {region: None}

    # Now try to run the CLONE_ASSET_APP num_retries + 1 times.
    curr_try = 0
    while curr_try <= num_retries:
        cmd = ['dx', 'run', CLONE_ASSET_APP_NAME, '--project', project, '-iurl=' + url, '-irecord_name=' + record_name]
        cmd += ['-iasset_file_name=' + asset_file_name, '-iasset_properties=' + json.dumps(asset_properties), '--brief']
        job = subprocess.check_output(cmd).strip()
        print('{0}: {1}'.format(region, job), file=sys.stderr)
        try:
            cmd = 'dx wait {0} '.format(job)
            subprocess.check_output(cmd, shell=True)
        except:
            pass
        cmd = 'dx describe {0} --json '.format(job)
        job_desc = json.loads(subprocess.check_output(cmd, shell=True).strip())

        if job_desc['state'] == 'done':
            q.put(region)
            return {region: job_desc['output']['asset_bundle']}

        curr_try += 1

    return {region: None}


def clone_asset(record_id, regions, num_retries=0, priority=None):
    """
    This function will attempt to clone the given record into all of the given regions.
    It will return a dictionary with the regions as keys and the record-ids of the
    corresponding asset as the values.  If an asset is not able to be created in a given
    region, the value will be set to None.
    """
    # Get the asset record
    record = dxpy.DXRecord(record_id)
    fid = record.get_details()['archiveFileId']['$dnanexus_link']
    curr_region = dxpy.describe(record.project)['region']

    # Only run once per region
    regions = set(regions) - set([curr_region])
    app_supported_regions = set(CLONE_ASSET_APP.describe()['regionalOptions'].keys())
    if len(regions - app_supported_regions) > 0:
        print('Currently no support for the following region(s): [{0}]'.format(', '.join(regions - app_supported_regions)), file=sys.stderr)
        sys.exit(1)

    # Get information about the asset
    record_name = record.name
    asset_properties = record.get_properties()
    asset_properties['cloned_from'] = record_id
    asset_file_name = dxpy.describe(fid)['name']
    url = dxpy.DXFile(fid).get_download_url(preauthenticated=True,
                                            project=dxpy.DXFile.NO_PROJECT_HINT,
                                            duration=URL_DURATION)[0]

    # Fire off a clone process for each region
    pool = multiprocessing.Pool(len(regions))
    manager = multiprocessing.Manager()
    q = manager.Queue()
    clone_asset_func = functools.partial(_clone_asset_into_region,
                                         record_name=record_name, q=q,
                                         asset_properties=asset_properties,
                                         asset_file_name=asset_file_name,
                                         url=url, num_retries=num_retries,
                                         priority=priority)
    results = pool.map_async(clone_asset_func, regions)

    # Get and return the results
    remaining_regions = regions
    print('Waiting on region(s): {0} '.format(' '.join(remaining_regions)))
    while True:
        if results.ready():
            break
        else:
            if q.qsize() > 0:
                for i in xrange(q.qsize()):
                    received = set([q.get()])
                    remaining_regions = remaining_regions - received
                print('\nWaiting on region(s): {0} '.format(' '.join(remaining_regions)))
            sys.stdout.write('.')
            sys.stdout.flush()
            time.sleep(SLEEP_TIME)

    print('\nDone')
    results = results.get()
    record_ids = {}
    for result in results:
        for region in result:
            if result[region] is None:
                record_ids[region] = None
            else:
                record_ids[region] = result[region]['$dnanexus_link']

    return record_ids


def main(record, regions, num_retries=0, priority=None):
    record_ids = clone_asset(record, regions, num_retries, priority)

    for region in record_ids:
        record_id = 'Failed' if record_ids[region] is None else record_ids[region]
        print('{0}:\t{1}'.format(region, record_id))


if __name__ == '__main__':
    args = _parse_args()
    main(args.record, args.regions, args.num_retries, args.priority)
