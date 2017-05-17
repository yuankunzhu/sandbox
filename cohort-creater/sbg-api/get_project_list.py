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


def get_project_list():
    for item in api.projects.query().all():
        list.append({'id': item.id, 'name': item.name})
    return list

print list
