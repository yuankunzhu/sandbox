#!/usr/bin/python
# -------------------------------------- #
# facade for sbg api on cohort-some
# -------------------------------------- #


# modules ====================================================================
import ConfigParser, requests, warlock, json
import sevenbridges as sbg
from flask import *


# configs =====================================================================
configs = ConfigParser.ConfigParser()
configs.read('config.ini')
url = configs.get('api', 'url')


# schema =====================================================================
payload = warlock.model_factory({
    "name": "payload",
    "additionalProperties": False,
    "properties": {
        "token": {"type": "string"},
        "username": {"type": "string"},
        "project": {"type": "string"},
        "projects": {"type": "array", "items": {"type": "string"}},
        "files": {"type": "array", "items": {"type": "string"}},
        "pipeline": {"type": "string"}
    }})

# input = payload(token='', projects=[], files=[], pipeline='')
payload_input = payload()
app = Flask(__name__)


# functions ==================================================================
# cohortsome page
@app.route('/cohortsome')
def cohortsome():
    return render_template('cohortsome.html')


# error handler
@app.errorhandler(500)
def special_exception_handler(error):
    return jsonify(codes=500)


# get user info/valid user token
@app.route('/api/v0/users')  # list user info
def get_user_info():
    payload_input.token = request.args.get('token')
    response = requests.request(
        method='GET',
        url=url + '/user',
        headers={"X-SBG-Auth-Token": payload_input.token})
    response_json = json.loads(response.content)
    return jsonify(response=response_json, codes=response.status_code)


# list avaibale projects
@app.route('/api/v0/users/sbg-project')
def list_sbg_project():
    list = []
    api = sbg.Api(url=url, token=request.args.get('token'))
    for item in api.projects.query().all():
        list.append({'id': item.id, 'name': item.name})
        for member in item.get_members().all():
            if member.username == api.users.me().username:
                member.permissions
    return jsonify(response=list)


# list files in files
@app.route('/api/v0/users/sbg-project/files')
def list_sbg_file():
    list = []
    api = sbg.Api(url=url, token=request.args.get('token'))
    for item in api.files.query(project=request.args.get('project')).all():
        list.append({'id': item.id, 'name': item.name})
    return jsonify(response=list)


# run task
@app.route('/api/v0/users/sbg-project/tasks')
def run_task():
    list = []

app.run(debug=True)
# ============================================================================
