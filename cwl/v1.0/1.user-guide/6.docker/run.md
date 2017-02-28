```
cwl-runner docker.cwl docker-job.yml
/usr/local/bin/cwl-runner 1.0.20170224141733
Resolved 'docker.cwl' to 'file:///Users/zhuy/Documents/chopwork/github/sandbox/cwl/v1.0/1.user-guide/docker/docker.cwl'
['docker', 'pull', 'node:slim']
slim: Pulling from library/node
5040bd298390: Pull complete
fce5728aad85: Pull complete
c8b82ee27706: Pull complete
2d64b1895c86: Pull complete
c3656194a8a8: Pull complete
Digest: sha256:9cbad37d57b2186a3a4a489944119ad42dc5424e09d0188dea3f85ccc94f0500
Status: Downloaded newer image for node:slim
[job docker.cwl] /var/folders/z8/pqc4mk8x1nzd8cy74zrlv67dx2txhk/T/tmpGjWA7o$ docker \
    run \
    -i \
    --volume=/Users/zhuy/Documents/chopwork/github/sandbox/cwl/v1.0/1.user-guide/docker/hello.js:/private/var/lib/cwl/stgddcf593c-be31-47f7-8381-dea36c10ca96/hello.js:ro \
    --volume=/private/var/folders/z8/pqc4mk8x1nzd8cy74zrlv67dx2txhk/T/tmpGjWA7o:/private/var/spool/cwl:rw \
    --volume=/private/var/folders/z8/pqc4mk8x1nzd8cy74zrlv67dx2txhk/T/tmpDnO97o:/tmp:rw \
    --workdir=/private/var/spool/cwl \
    --read-only=true \
    --user=2049799698 \
    --rm \
    --env=TMPDIR=/tmp \
    --env=HOME=/private/var/spool/cwl \
    node:slim \
    node \
    /private/var/lib/cwl/stgddcf593c-be31-47f7-8381-dea36c10ca96/hello.js
Hello World
[job docker.cwl] completed success
{}
Final process status is success
```