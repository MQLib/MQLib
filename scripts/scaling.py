# A script to test the memory scaling properties of the algorithms. This script
# creates Erdos-Renyi graphs with weights randomly selected from the set {-1, 1}.
# One set of graphs is dense and the other is sparse. It tests all the heuristics
# on graphs of specified sizes, outputting the memory usage.
import csv
import os.path
import platform
import random
import subprocess
import sys

# Create a complete graph with random edge weights selected from {-1, 1}
def createCompleteGraph(size):
    random.seed(144)
    with open('complete.txt', 'w') as fp:
        fp.write(str(size) + " " + str(size*(size-1)/2) + "\n")
        for x in range(size-1):
            for y in range(x+1, size):
                fp.write(str(x+1) + " " + str(y+1) + " " + str(2*random.randint(0, 1)-1) + "\n")

# Create an Erdos-Renyi graph with random edge weights selected from {-1, 1}
# and edge probabilities that yield an expected degree of 5 for each node.
def createERGraph(size):
    random.seed(144)
    p = 5.0 / (size-1)
    edges = []
    for x in range(size-1):
        for y in range(x+1, size):
            if random.random() < p:
                edges.append((x, y))
    with open('ER.txt', 'w') as fp:
        fp.write(str(size) + " " + str(len(edges)) + "\n")
        for x, y in edges:
            fp.write(str(x+1) + " " + str(y+1) + " " + str(2*random.randint(0, 1)-1) + "\n")

# Determine the runtime limit for the indicated graph
def getRuntimeLimit(graphName):
    torun = ["timeout", "10000", "../bin/MQLib", "-fM",
             graphName, "-h", "BASELINE", "-r", "1500", "-s", "144"]
    p = subprocess.Popen(torun, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    baseline_output = p.stdout.read().decode("utf-8")
    return baseline_output.split(",")[4]


##############################
##############################
minERGraphSize = 1000
maxCompleteGraphSize = 3000
if len(sys.argv) < 3 or not all([x.isdigit() for x in sys.argv[2:]]):
    print('Usage: python scaling.py heuristics.txt n1 n2 n3 ...')
    exit(0)

if not os.path.exists("../bin"):
    print("scaling.py must be run from the scripts folder")
    exit(0)
if not os.path.exists("../bin/MQLib"):
    print("You need to run `make` in the main folder before running scaling.py")
    exit(0)

# Load the heuristics to evaluate
with open(sys.argv[1], 'r') as fp:
    heuristics = [x.strip() for x in fp]

writer = csv.writer(sys.stdout)
writer.writerow(['heuristic', 'graph', 'size', 'runtime', 'memusg'])

for size in [int(x) for x in sys.argv[2:]]:
    runs = []
    if size <= maxCompleteGraphSize:
        createCompleteGraph(size)
        runs.append(('complete.txt', getRuntimeLimit('complete.txt')))
    if size >= minERGraphSize:
        createERGraph(size)
        runs.append(('ER.txt', getRuntimeLimit('ER.txt')))

    for graphName, runtime in runs:
        for heuristic in heuristics:
            if platform.system() == 'Darwin':
                torun = ['/usr/bin/time', '-l', '../bin/MQLib', '-fM',
                         graphName, '-h', heuristic, '-r', runtime, '-s', '144',
                         '-nv']
                p = subprocess.Popen(torun, stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT)
                baseline_output = p.stdout.read().decode("utf-8")
                for q in baseline_output.split("\n"):
                    if q.find("maximum resident set size") >= 0:
                        writer.writerow([heuristic, graphName, size, runtime,
                                         q.strip().split()[0]])
            else:
                torun = ['/usr/bin/time', '-v', '../bin/MQLib', '-fM',
                         graphName, '-h', heuristic, '-r', runtime, '-s', '144',
                         '-nv']
                p = subprocess.Popen(torun, stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT)
                baseline_output = p.stdout.read().decode("utf-8")
                for q in baseline_output.split("\n"):
                    if q.find("Maximum resident set size (kbytes):") >= 0:
                        writer.writerow([heuristic, graphName, size, runtime,
                                         q.strip().split()[5]])
