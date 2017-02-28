Normally, input files are located in a read-only directory separate from the output directory. This causes problems if the underlying tool expects to write its output files alongside the input file in the same directory. You use `InitialWorkDirRequirement` to stage input files into the output directory. In this example, we use a Javascript expression to extract the base name of the input file from its leading directory path.
```
/usr/local/bin/cwl-runner 1.0.20170224141733
Resolved 'linkfile.cwl' to 'file:///Users/zhuy/Documents/chopwork/github/sandbox/cwl/v1.0/1.user-guide/14.staging-inputs/linkfile.cwl'
['docker', 'pull', 'java:7']
7: Pulling from library/java
Digest: sha256:c0b61b62639124aa838dc755c5a9d57c072f762b71b170281927399a14db4652
Status: Downloaded newer image for java:7
[job linkfile.cwl] /var/folders/z8/pqc4mk8x1nzd8cy74zrlv67dx2txhk/T/tmpsYt7Y8$ docker \
    run \
    -i \
    --volume=/Users/zhuy/Documents/chopwork/github/sandbox/cwl/v1.0/1.user-guide/14.staging-inputs/Hello.java:/private/var/lib/cwl/stg10bc3e2c-a510-45a3-b965-17651d9a522f/Hello.java:ro \
    --volume=/private/var/folders/z8/pqc4mk8x1nzd8cy74zrlv67dx2txhk/T/tmpsYt7Y8:/private/var/spool/cwl:rw \
    --volume=/private/var/folders/z8/pqc4mk8x1nzd8cy74zrlv67dx2txhk/T/tmpgQhEy_:/tmp:rw \
    --workdir=/private/var/spool/cwl \
    --read-only=true \
    --user=2049799698 \
    --rm \
    --env=TMPDIR=/tmp \
    --env=HOME=/private/var/spool/cwl \
    java:7 \
    javac \
    Hello.java
[job linkfile.cwl] completed success
{
    "classfile": {
        "checksum": "sha1$e68df795c0686e9aa1a1195536bd900f5f417b18",
        "basename": "Hello.class",
        "location": "file:///Users/zhuy/Documents/chopwork/github/sandbox/cwl/v1.0/1.user-guide/14.staging-inputs/Hello.class",
        "path": "/Users/zhuy/Documents/chopwork/github/sandbox/cwl/v1.0/1.user-guide/14.staging-inputs/Hello.class",
        "class": "File",
        "size": 184
    }
}
Final process status is success
```