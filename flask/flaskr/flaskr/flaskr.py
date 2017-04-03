# all the imports
import os
import sqlite3
from flask import Flask
# import request, session, g, redirect, url_for, abort, render_template, flash

local_cfg = dict(
    USERNAME='zhuyk'
)

# app = Flask(__name__)
# app.config.update(local_cfg)
# print app.config

app = Flask(__name__)  # create the application instance :)
app.config.from_object(__name__)  # load config from this file , flaskr.py
# Load default config and override config from an environment variable


app.config.update(dict(
    DATABASE=os.path.join(app.root_path, 'flaskr.db'),
    SECRET_KEY='development key',
    USERNAME='admin',
    PASSWORD='default'
))
app.config.from_envvar('FLASKR_SETTINGS', silent=True)
# app.config.from_pyfile


def connect_db():
    """Connects to the specific database."""
    rv = sqlite3.connect(app.config['DATABASE'])
    rv.row_factory = sqlite3.Row
    return rv
