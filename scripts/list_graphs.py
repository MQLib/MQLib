import boto
import os
import sys

conn = boto.connect_s3(anon=True)
b = conn.get_bucket("mqlibinstances", validate=False)

for key in b.list():
    print(key.name)
