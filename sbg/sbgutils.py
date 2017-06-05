import logging
import json
import requests
import os
import pandas as pd
import sqlalchemy as sa

def get_sbg_metadata(dataset, tokenfile, cavatica_base, metadata_fields, outputdir='.'):
    filename = os.path.join(outputdir, '{}.sbg_metadata.csv'.format(dataset.split('/')[-1]))
    json_filename = os.path.join(outputdir, '{}.sbg_metadata.json'.format(dataset.split('/')[-1]))
    if os.path.isfile(filename):
        sbg_df = pd.read_csv(filename, dtype={'aliquot_id': str, 'file_name': str, 'sample_id': str})
    else:
        sbg_df, jsons = get_files_aliquots_sbg(dataset, tokenfile, cavatica_base, metadata_fields)
        sbg_df.to_csv(filename, index=False)
        with open(json_filename, 'w') as json_file:
            json_file.write(json.dumps(jsons, indent=2))

    #create new row for each "pooled" aliquot
    logging.info("pre pooling sbg_df.shape: {}".format(sbg_df.shape))
    sbg_noaliquot_df = sbg_df[pd.isnull(sbg_df['aliquot_id'])]

    sbg_df = sbg_df[pd.notnull(sbg_df['aliquot_id'])]
    sbg_df['aliquot_ids'] = sbg_df['aliquot_id'].apply(lambda x: x.split('-'))

    s = sbg_df.apply(lambda x: pd.Series(x['aliquot_ids']),axis=1).stack().reset_index(level=1, drop=True)
    s.name = 'aliquot_id'
    sbg_df = sbg_df.drop('aliquot_id', axis=1).join(s)
    logging.info("after pooling sbg_df.shape: {}".format(sbg_df.shape))

    return sbg_df, sbg_noaliquot_df

def get_files_aliquots_sbg(dataset, tokenfile, cavatica_base, metadata_fields, limit=100):
    token = open(tokenfile).readline().strip()
    headers = {'X-SBG-Auth-Token': token, 'Content-type': 'application/json'}
    ds_file_url = '{}/{}?limit={}&dataset={}'.format(cavatica_base, 'files', limit, dataset)
    items = {'sbg_file_url': [], 'file_size': [], 'file_name': [], 'dataset' : dataset.split('/')[-1]}
    for field in metadata_fields:
        items[field] = []

    jsons = []
    while ds_file_url:
        logging.info("ds_file_url: {}".format(ds_file_url))

        ds_resp = requests.get(ds_file_url, headers=headers)
        files = ds_resp.json()

        logging.debug("resp: {}".format(files))
        for f in files['items']:
            details_resp = requests.get(f['href'], headers=headers)
            file_details = details_resp.json()
            jsons.append(file_details)
            logging.debug(json.dumps(file_details, indent=2))
            for field in metadata_fields:
                if field not in file_details['metadata']:
                    items[field].append(None)
                else:
                    items[field].append(file_details['metadata'][field])

            items['sbg_file_url'].append(f['href'])
            items['file_size'].append(file_details['size'])
            items['file_name'].append(f['name'])

        next_url = None
        for link in files['links']:
            if link['rel'] == 'next':
                next_url = link['href']

        ds_file_url = next_url

    sbg_df = pd.DataFrame(items)
    return sbg_df, jsons

class HarvestDB:

    def __init__(self, host, database, username, password):
        self.host = host
        self.database = database
        self.username = username
        self.password = password

        if username:
            self.engine = sa.create_engine('postgresql://{}:{}@{}/{}'.format(username, password, host, database))
        else:
            self.engine = sa.create_engine('postgresql://{}/{}'.format(host, database))

        self.conn = self.engine.connect()
        self.meta = sa.MetaData()


    def get_table_as_df(self, tablename):
        table = sa.Table(tablename, self.meta,
                         autoload=True, autoload_with=self.engine)
        s = sa.sql.select([table])
        result = self.conn.execute(s)
        df = pd.DataFrame(result.fetchall())
        df.columns = result.keys()
        return df