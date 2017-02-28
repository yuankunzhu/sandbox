```bash
$ rm hello.tar || true && touch goodbye.txt && tar -cvf hello.tar goodbye.txt
$ cwl-runner tar-param.cwl tar-param-job.yml
[job 139868145165200] $ tar xf /home/example/hello.tar goodbye.txt
Final process status is success
{
"example_out": {
  "location": "goodbye.txt",
  "size": 24,
  "class": "File",
  "checksum": "sha1$dd0a4c4c49ba43004d6611771972b6cf969c1c01"
  }
}
```