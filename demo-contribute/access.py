import os, csv, time
import tempfile, shutil
import pandas as pd
import sevenbridges as sbg

ctrl_projects = [
    'cavatica/childrens-cancer-atlas',
    'cavatica/pacific-pediatric-neuro-oncology-consortium',
    'cavatica/pediatric-cardiac-genomics-consortium',
    'cavatica/childhood-brain-tumor-tissue-consortium',
    'cavatica/john-maris-laboratory',
    'cavatica/pediatric-cancer-dream-team',
    'cavatica/22q11-deletion-syndrome-project'
]

token = os.environ['CAVATICA_TOKEN']
url = 'https://cavatica-api.sbgenomics.com/v2/'

def pending_user_dict(type, member, df, datasetName):
    try:
        user = df.loc[member.username]
        fullname = str(user['first name']) + ' ' + str(user['last name'])
        affiliation = user.affiliation
        email = user.email
    except:
        affiliation = '?'
        email = '?'
        fullname = '?'
    return ({
        'DataCollectionType': type,
        'DataCollectionName': datasetName,
        'Username': member.username,
        'Fullname': fullname,
        'Affiliation': affiliation,
        'Email': email,
        'WriteAccess': member.permissions['write'],
        'CopyAccess': member.permissions['copy'],
        'ExecutAccess': member.permissions['execute'],
        'AdminAccess': member.permissions['admin']})


def main():

    api = sbg.Api(url=url, token=token)
    stats_files = api.files.query(project='marko_sbg/stats').all()
    user_stats = [file for file in stats_files if 'users' in file.name]

    tmpdir = tempfile.mkdtemp()
    df_list = []

    for each_stats in user_stats:
        stats_csv = os.path.join(tmpdir, each_stats.name)
        each_stats.download(path=stats_csv)
        df_list.append(pd.read_csv(stats_csv))

    df = pd.concat(df_list).drop_duplicates(subset='username', keep='last')
    df.set_index('username')
    shutil.rmtree(tmpdir)

    pending_user = []

    for dataset in api.datasets.query():
        for member in dataset.members().all():
            pending_user.append(
                pending_user_dict(
                    'DATASETS', member, df, dataset.name))

    for ctrl_project in ctrl_projects:
        project = api.projects.get(id=ctrl_project)
        members = project.get_members().all()
        for member in members:
            pending_user.append(
                pending_user_dict(
                    'CONTROLLED', member, df, dataset.name))

    keys = pending_user[0].keys()
    timestamp = time.strftime("%Y-%m-%d.%H%M%S")
    output_csv = 'pending_user.' + timestamp + '.csv'

    with open(output_csv, 'wb') as fo:
        dict_writer = csv.DictWriter(fo, keys)
        dict_writer.writeheader()
        dict_writer.writerows(pending_user)


if __name__ == '__main__':
    main()
