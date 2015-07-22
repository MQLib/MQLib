import os.path
import subprocess
import sys

if len(sys.argv) != 2:
    print "Usage: python INSTALL.py http://link/to/git/repo.git"
    exit(1)

# Run installation script until it works (sometimes it takes more than once due to
# apt-get servers not working, etc.)
iteration = 0
while True:
    iteration += 1
    try:
        p = subprocess.Popen(["timeout", "360", "./INSTALL.sh", str(iteration),
                              sys.argv[1]], stdout=None, stderr=None)
        p.wait()
    except:
        with open("failure.txt", "w") as f:
            f.write("INSTALL.py exited in failure\n")
        exit(1)

    if os.path.isfile("MQLib/bin/MQLib"):
        break
