import argparse
import json
import requests
import time
'''
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
    '''
def pdb_to_dssp(pdb_file_path, rest_url):
    # Read the pdb file data into a variable
    files = {'file_': open(pdb_file_path, 'rb')}

    url_create = '{}api/create/pdb_file/dssp/'.format(rest_url)
    r = requests.post(url_create, files=files)
    r.raise_for_status()

    job_id = json.loads(r.text)['id']
    print ("Job submitted successfully. Id is: '{}'".format(job_id))

    # Loop until the job running on the server has finished, either successfully
    # or due to an error.
    ready = False
    while not ready:
        # Check the status of the running job. If an error occurs an exception
        # is raised and the program exits. If the request is successful, the
        # status is returned.
        url_status = '{}api/status/pdb_file/dssp/{}/'.format(rest_url,
                                                                  job_id)
        r = requests.get(url_status)
        r.raise_for_status()

        status = json.loads(r.text)['status']
        print ("Job status is: '{}'".format(status))

        # If the status equals SUCCESS, exit out of the loop by changing the
        # condition ready. This causes the code to drop into the `else` block
        # below.
        #
        # If the status equals either FAILURE or REVOKED, an exception is raised
        # containing the error message. The program exits.
        #
        # Otherwise, wait for five seconds and start at the beginning of the
        # loop again.
        if status == 'SUCCESS':
            ready = True
        elif status in ['FAILURE', 'REVOKED']:
            raise Exception(json.loads(r.text)['message'])
        else:
            time.sleep(5)
    else:
        # Requests the result of the job. If an error occurs an exception is
        # raised and the program exits. If the request is successful, the result
        # is returned.
        url_result = '{}api/result/pdb_file/dssp/{}/'.format(rest_url,
                                                                  job_id)
        r = requests.get(url_result)
        r.raise_for_status()
        result = json.loads(r.text)['result']

        # Return the result to the caller, which prints it to the screen.
        return result
