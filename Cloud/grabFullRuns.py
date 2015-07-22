# python grabFullRuns.py
import csv
import os
import sys
import CloudSetup

if len(sys.argv) != 2:
    print "Usage: python grabFullRuns.py runtimes.csv"
    exit(1)

runtimes = {}  # graphname -> runtime limit
with open(sys.argv[1], "rU") as f:
    r = csv.reader(f)
    if r.next() != ["graphname", "runtime"]:
        print "Illegal header for", sys.argv[1]
        exit(1)
    for line in r:
        runtimes[line[0]] = line[1]

suffixes = [""] + [str(x) for x in range(1, 10000)]
resultsFile = None
errorsFile = None
for s in suffixes:
    if not os.path.exists("../data/results" + s + ".csv") and not os.path.exists("../data/errors" + s + ".txt"):
        resultsFile = "../data/results" + s + ".csv"
        errorsFile = "../data/errors" + s + ".txt"
        break
if resultsFile is None or errorsFile is None:
    print "Could not allocate results or errors file"
    exit(1)
print "Outputting results to", resultsFile
print "Outputting errors to", errorsFile

sdb, dom = CloudSetup.setup_sdb_domain("mqlib-domain")
rs = dom.select('select * from `mqlib-domain`')
dat = {}  # graphname -> output
writer = csv.writer(open(resultsFile, "w"))
writer.writerow(["timestamp", "graphname", "heuristic", "limit", "objective",
                 "runtime"])
errorOut = open(errorsFile, "w")
for result in rs:
    graphname = result["graphname"].strip()
    heuristic = result["heuristic"].strip()
    timestamp = result["timestamp"].strip()
    output = ""
    for idx in range(1000):
        key = "output" + str(idx).zfill(3)
        if key in result:
            output += result[key]
        else:
            break  # We have exhausted all the output fields
    output = output.strip()

    # The heuristic history is between the []'s and is of the form
    # objective1:runtime1;objective2:runtime2;...
    start = output.find("[")
    end = output.find("]")
    if start < 0 or end < 0 or start >= end-1:
        # Error in output; just report the solution of 0 at 0 seconds and log
        # to the errors file.
        if not graphname in runtimes:
            rtlim = "-1"  # No runtime limit in the runtimes files...
        else:
            rtlim = runtimes[graphname]
        writer.writerow([timestamp, graphname, heuristic, rtlim, "0", "0"])
        errorOut.write("************ Oddly formatted heuristic output for " +
                       graphname + " (" + heuristic + ")\n")
        errorOut.write("timestamp: " + timestamp + "\n")
        errorOut.write(output + "\n")
    else:
        limit = output.split(",")[0]
        part = output[(start+1):end]
        for piece in part.split(";"):
            obj = piece.split(":")[0]
            rt = piece.split(":")[1]
            if float(obj) == 0:
                rt = "0"  # We know solution 0 from the very beginning
            if float(rt) > float(limit):
                continue  # We don't process results after RT limit
            writer.writerow([timestamp, graphname, heuristic, limit, obj, rt])
