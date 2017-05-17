#!/usr/bin/python
# -------------------------------------- #
# sbg api for cohort-run to cbio
# -------------------------------------- #

# modules ====================================================================
import ConfigParser
import sys, json
import sevenbridges as sbg

# configs =====================================================================
configs = ConfigParser.ConfigParser()
configs.read('config.ini')
url = configs.get('api', 'url')

# intial api & payload ========================================================
with open(sys.argv[1]) as data_file:
    data = json.load(data_file)
api = sbg.Api(url=url, token=data['sbg_token'])


# functions ===================================================================
# get the project list with permissions
def get_project_list():
    project_list = []
    for item in api.projects.query().all():
        for member in item.get_members().all():
            if member.username == api.users.me().username:
                project_list.append({
                    'id': item.id,
                    'name': item.name,
                    'exe': member.permissions['execute']})
    return project_list


# get the file list with meta
def get_file_list():
    file_list = []
    for item in api.files.query(project=data['project_id']).all():
        if item.metadata is None:
            meta = False
        elif(item.metadata['experimental_strategy'] is None):
            meta = False
            seqtype = False
        elif(item.metadata['paired_end'] and item.metadata['sample_id'] is not None):
            meta = True
            seqtype = item.metadata['experimental_strategy']
        file_list.append({
            'id': item.id,
            'name': item.name,
            'metadata': meta,
            'seqtype': seqtype})
    return file_list


# cohort-run
def cohort_cbio_run():
    files = data['files']
    inputs = {}
    inputs['readFilesIn'] = api.files.query(project=data['project_id'], names=files)
    # inputs['readFilesIn'] = api.files.query(project='yuankun/helloworld', metadata={'sample_id': 'ERR315451'})
    inputs['indextar'] = api.files.get(id='5800dc01e4b0144f66f20f52').copy(project=data['project_id'])
    inputs['genome'] = api.files.get(id='5800dc01e4b0144f66f20f44').copy(project=data['project_id'])

    inputs['outFileNamePrefix'] = ''
    inputs['sample_name'] = ''

    task = api.tasks.create(
        name=data['cohort_name'],
        project=data['project_id'],
        app=data['pipeline_id'],
        inputs=inputs)


# main ========================================================================
if __name__ == '__main__':
    # get_project_list()
    # get_file_list()
    cohort_cbio_run()
