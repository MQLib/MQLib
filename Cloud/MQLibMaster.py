import subprocess
import sys

# Validate command-line arguments
if len(sys.argv) < 2 or (not (sys.argv[1] == "METRICS" and len(sys.argv) == 3) and not (sys.argv[1] == "FULL" and len(sys.argv) == 7 and sys.argv[3].isdigit() and all([x.isdigit() for x in sys.argv[4].split("_")]) and sys.argv[5].lstrip("-").isdigit() and sys.argv[6].lstrip("-").isdigit())):
    print("Usage:\n  python MQLibMaster.py METRICS tag\n    [[or]]\n  python MQLibMaster.py FULL tag #ITERFORBASELINE SEEDS_SEPARATED_BY_UNDERSCORES MINSECONDS MAXSECONDS")
    exit(1)

# Run until it tells us that we're done
while True:
    if sys.argv[1] == "METRICS":
        p = subprocess.Popen(["python", "MQLibRunner.py", sys.argv[1], sys.argv[2]],
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    else:
        p = subprocess.Popen(["python", "MQLibRunner.py", sys.argv[1], sys.argv[2], sys.argv[3],
                              sys.argv[4], sys.argv[5], sys.argv[6]],
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    for line in p.stdout:
        sys.stdout.write(line)
    p.wait()

    # MQLibRunner.py will terminate this EC2 node if it completes successfully,
    # so if we're still running then it must have failed. We'll just kick
    # it off again at the top of the loop.
