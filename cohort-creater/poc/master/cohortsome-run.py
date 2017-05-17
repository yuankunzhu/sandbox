#!/usr/bin/python
# -------------------------------------- #
# facade for sbg api on cohort-some
# -------------------------------------- #


# modules ====================================================================
import ConfigParser, requests, warlock, json
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
        "token": {
            "type": "string"
        },
        "username": {
            "type": "string"
        },
        "project": {
            "type": "string"
        },
        "projects": {
            "type": "array",
            "items": {"type": "string"}
        },
        "files": {
            "type": "array",
            "items": {"type": "string"}
        },
        "pipeline": {"type": "string"}
    }})

# input = payload(token='', projects=[], files=[], pipeline='')
payload_input = payload()
app = Flask(__name__)


# functions ==================================================================
@app.route('/cohortsome')  # list user info
def cohortsome():
    return render_template('cohortsome.html')


# functions ==================================================================
@app.errorhandler(500)
def special_exception_handler(error):
    return jsonify(codes=500)


# functions ==================================================================
@app.route('/api/v0/users')  # list user info
def get_user_info():
    payload_input.token = request.args.get('token')
    response = requests.request(
        method='GET',
        url=url + '/user',
        headers={"X-SBG-Auth-Token": payload_input.token}
    )
    response_json = json.loads(response.content)
    return jsonify(response=response_json, codes=response.status_code)


# functions ==================================================================
@app.route('/api/v0/users/sbg-project')  # list avaibale projects
def list_sbg_project():
    payload_input.token = request.args.get('token')
    response = requests.request(
        method='GET',
        url=url + '/projects',
        headers={"X-SBG-Auth-Token": payload_input.token}
    )
    response_json = json.loads(response.content)
    return jsonify(response=response_json, codes=response.status_code)


# functions ==================================================================
@app.route('/api/v0/users/sbg-project/files')  # list files in project
def list_sbg_file():
    payload_input.token = request.args.get('token')
    payload_input.username = request.args.get('user')
    payload_input.project = request.args.get('project')
    response = requests.request(
        method='GET',
        url=url + '/projects/' + payload_input.project + '/files',
        headers={"X-SBG-Auth-Token": payload_input.token}
    )
    response_json = json.loads(response.content)
    return jsonify(response=response_json, codes=response.status_code)


app.run(debug=True)
# ============================================================================
