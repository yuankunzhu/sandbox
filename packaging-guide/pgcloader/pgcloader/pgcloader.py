import sys, json, csv
import requests, time
import ConfigParser

config = sys.argv[1]

configs = ConfigParser.ConfigParser()
configs.read(config)

token = configs.get('default', 'token')
url = configs.get('default', 'url')
manifest = configs.get('parameter', 'manifest')
project = configs.get('parameter', 'project')
sleeptime = configs.get('parameter', 'sleep')

fields = list(csv.reader(open(manifest), delimiter='\t'))


def api(path, method='GET', data=None):
	data = json.dumps(data) if isinstance(data, dict) else None
	headers = {"X-SBG-Auth-Token": token, "Content-type": "application/json"}
	response = requests.request(method, path, data=data, headers=headers)
	response_json = json.loads(response.content) if response.content else {}
	return response_json


def dump(content):
	print json.dumps(content, indent=4)


def main():
	rate = api(path=url + 'rate_limit')
	files = []
	metadata = {}
	for field in fields:
		field_split = field[0].split('/')
		bucket = field_split[2]
		volume = configs.get('bucket', bucket)
		location = '/'.join(field_split[3:])
		body_import = {
			"source": {"volume": volume, "location": location},
			"destination": {"project": project},
			"overwrite": True
		}
		data_import = api(
			path=url + 'storage/imports',
			method='POST',
			data=body_import
		)
		files.append(data_import['id'])
		metadata[data_import['id']] = field[1]
		rate['rate']['remaining'] -= 1
		while rate['rate']['remaining'] <= 1 and rate['rate']['reset'] > int(time.time()):
			time.sleep(5)
		if rate['rate']['reset'] < int(time.time()):
			rate = api(path=url + 'rate_limit')
	print "all file imported; in " + sleeptime + " sec, will start patching metadata..."
	time.sleep(int(sleeptime))
	for file in files:
		import_detail = api(path=url + 'storage/imports/' + file)
		meta_body = metadata[file]
		meta_url = url + 'files/' + import_detail['result']['id'] + '/metadata'
		api(path=meta_url, method='PUT', data=json.loads(meta_body))
		rate['rate']['remaining'] -= 1
		while rate['rate']['remaining'] <= 1 and rate['rate']['reset'] > int(time.time()):
			time.sleep(5)
		if rate['rate']['reset'] < int(time.time()):
			rate = api(path=url + 'rate_limit')

if __name__ == '__main__':
	main()
	print "all done"
