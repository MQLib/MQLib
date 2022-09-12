import boto
import csv
import sys
import io
import zipfile

writer = csv.writer(sys.stdout)
writer.writerow(['fname', 'n', 'm', 'comments'])

conn = boto.connect_s3(anon=True)
b = conn.get_bucket("mqlibinstances", validate=False)
for k in b.list():
    name = k.key
    sys.stderr.write(name + " ...\n")
    contents = None
    for count in range(10):
        try:
            stringContents = k.get_contents_as_string()
            contents = zipfile.ZipFile(io.BytesIO(stringContents), "r")
            break
        except:
            sys.stderr.write("  Download iteration failed on " + name + "\n")
    if contents is None:
        sys.stderr.write("Failed to download " + name + "\n")
        exit(1)
    fnames = [x.filename for x in contents.infolist()]
    if len(fnames) != 1:
        sys.stderr.write("Bad zip file: " + name + "\n")
        exit(1)
    with contents.open(fnames[0]) as fp:
        first = None
        comments = []
        for line in fp:
            line = line.decode("utf-8")
            if line.strip().find("#") == 0:
                comments.append(line[1:].strip())
            elif first is None:
                first = line.strip()
        writer.writerow([name, first.split()[0], first.split()[1],
                         ' || '.join(comments)])

