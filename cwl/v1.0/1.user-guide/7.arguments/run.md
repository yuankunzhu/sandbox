```
$ cwl-runner arguments.cwl arguments-job.yml
/usr/local/bin/cwl-runner 1.0.20170224141733
Resolved 'arguments.cwl' to 'file:///Users/zhuy/Documents/chopwork/github/sandbox/cwl/v1.0/1.user-guide/7.arguments/arguments.cwl'
['docker', 'pull', 'java:7-jdk']
7-jdk: Pulling from library/java
5040bd298390: Already exists
fce5728aad85: Already exists
76610ec20bf5: Pull complete
60170fec2151: Pull complete
66b144b1d5b0: Pull complete
6263baad4f89: Pull complete
Digest: sha256:c0b61b62639124aa838dc755c5a9d57c072f762b71b170281927399a14db4652
Status: Downloaded newer image for java:7-jdk
[job arguments.cwl] /var/folders/z8/pqc4mk8x1nzd8cy74zrlv67dx2txhk/T/tmpF0Inqm$ docker \
    run \
    -i \
    --volume=/Users/zhuy/Documents/chopwork/github/sandbox/cwl/v1.0/1.user-guide/7.arguments/Hello.java:/private/var/lib/cwl/stg931a2e07-6f6a-4116-844a-50cf87ac1e48/Hello.java:ro \
    --volume=/private/var/folders/z8/pqc4mk8x1nzd8cy74zrlv67dx2txhk/T/tmpF0Inqm:/private/var/spool/cwl:rw \
    --volume=/private/var/folders/z8/pqc4mk8x1nzd8cy74zrlv67dx2txhk/T/tmpmSS49u:/tmp:rw \
    --workdir=/private/var/spool/cwl \
    --read-only=true \
    --user=2049799698 \
    --rm \
    --env=TMPDIR=/tmp \
    --env=HOME=/private/var/spool/cwl \
    java:7-jdk \
    javac \
    -d \
    /private/var/spool/cwl \
    /private/var/lib/cwl/stg931a2e07-6f6a-4116-844a-50cf87ac1e48/Hello.java
[job arguments.cwl] completed success
{
    "classfile": {
        "checksum": "sha1$e68df795c0686e9aa1a1195536bd900f5f417b18",
        "basename": "Hello.class",
        "location": "file:///Users/zhuy/Documents/chopwork/github/sandbox/cwl/v1.0/1.user-guide/7.arguments/Hello.class",
        "path": "/Users/zhuy/Documents/chopwork/github/sandbox/cwl/v1.0/1.user-guide/7.arguments/Hello.class",
        "class": "File",
        "size": 184
    }
}
Final process status is success
```

Here we use the arguments field to add an additional argument to the command line that isn't tied to a specific input parameter.

arguments: ["-d", $(runtime.outdir)]