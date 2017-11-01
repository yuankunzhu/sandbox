import boto3
# session = boto3.Session(profile_name='saml')
# s3 = session.resource('s3')

# object = s3.Object('kids-first-seq-data','x01/bcm/batch_01/HKNY3CCXX-2.hgv.bam')
# print object.e_tag

# object = s3.Object('kids-first-seq-data','x01/bcm/batch_01/HLKMYCCXX-3.hgv.bam')
# print object.e_tag

# object = s3.Object('kids-first-seq-data','x01/bcm/batch_01/Manifest.txt')
# print object.e_tag

session = boto3.Session(profile_name='default')
s3 = session.resource('s3')

object = s3.Object('cbttc.seq.data','NantWorks/008f259f-c6a4-49a1-b404-efe16afb3a8a.bam')
print object.e_tag