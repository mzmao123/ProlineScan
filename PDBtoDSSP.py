import argparse
import json
import requests
import time
def pdb_to_dssp(pdb_file_path,rest_url):
    files = {'file_':open(pdb_file_path,'rb')}
    url_create = '{}api/create/pdb_file/dssp/'.format(rest_url)
    r = requests.post(url_create, files=files)
    r.raise_for_status()
    job_id = json.loads(r.text)['id']
    # unsure if the status part is necessary, only thing we need is the dssp file
    """
    ready = False
    while not ready:
        url_status = '{}api/status/pdb_file/dssp/{}/'.format(rest_url,
                                                                  job_id)
        r = requests.get(url_status)
        r.raise_for_status()

        status = json.loads(r.text)['status']

        if status == 'SUCCESS':
            ready = True
        elif status in ['FAILURE', 'REVOKED']:
            raise Exception(json.loads(r.text)['message'])
        else:
            time.sleep(5)
    else:
    """
    url_result = '{}api/result/pdb_file/dssp/{}/'.format(rest_url,job_id)
    r = requests.get(url_result)
    r.raise_for_status()
    result = json.loads(r.text)['result']

    return result
