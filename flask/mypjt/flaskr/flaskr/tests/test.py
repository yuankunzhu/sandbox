from flask import Flask
import os


# print __name__
app = Flask(__name__)
print os.path.join(app.root_path, 'flaskr.db')
