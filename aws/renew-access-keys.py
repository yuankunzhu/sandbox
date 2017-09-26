import sys, os, json
import ConfigParser
import boto3

config_file = os.path.join(os.path.expanduser('~'), '.aws/credentials')
configs = ConfigParser.ConfigParser()
configs.read(config_file)

username = os.environ['AWS_USERNAME']
access_key = configs.get('default', 'aws_access_key_id')
secret_key = configs.get('default', 'aws_secret_access_key')

client = boto3.client('iam')
resource = boto3.resource('iam')

client.list_access_keys(Username=username)


for key in user.access_keys.all():
    if key.access_key_id != access_key:
        
