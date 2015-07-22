import subprocess
import sys

# Validate command-line arguments
if len(sys.argv) != 3 or not sys.argv[1] in ["METRICS", "FULL"]:
    print "Usage: python MQLibMaster.py METRICS|FULL tag"
    exit(1)

# Run until it tells us that we're done
while True:
    p = subprocess.Popen(["python", "MQLibRunner.py", sys.argv[1], sys.argv[2]],
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    for line in p.stdout:
        sys.stdout.write(line)
    p.wait()

    # MQLibRunner.py will terminate this EC2 node if it completes successfully,
    # so if we're still running then it must have failed. We'll just kick
    # it off again at the top of the loop.
