import boto
import sys

if len(sys.argv) != 2:
    print("Usage: python downloadGraph.py graphName.txt")
    exit(1)

conn = boto.connect_s3(anon=True)
b = conn.get_bucket("mqlibinstances", validate=False)
k = boto.s3.key.Key(b)
k.key = sys.argv[1]
if sys.version_info > (3, 0):
    sys.stdout.buffer.write(k.get_contents_as_string())
else:
    sys.stdout.write(k.get_contents_as_string())
