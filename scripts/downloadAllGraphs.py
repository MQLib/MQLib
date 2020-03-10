import boto
import os
import sys

if len(sys.argv) != 2:
    print("Usage: python downloadAllGraphs.py outputFolder")
    exit(1)

conn = boto.connect_s3(anon=True)
b = conn.get_bucket("mqlibinstances", validate=False)
for k in b.list():
    name = k.key
    fname = sys.argv[1] + "/" + name
    if os.path.isfile(fname):
        print("[ Skipping ", name, "]")
        continue
    print(name, "...")
    success = False
    for count in range(10):
        try:
            k.get_contents_to_filename(fname)
            success = True
        except:
            print("  download failed")
        if success:
            break
    if not success:
        print("Failed to download", name)
        exit(1)
