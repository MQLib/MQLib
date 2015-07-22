import csv
import math
import subprocess
import CloudSetup

# Compute the runtime limit from the metrics.
def runtime(metrics):
    rt = 0.59 * math.exp(metrics["log_n"])
    if rt < 120:
        return 120
    elif rt > 1200:
        return 1200
    else:
        return rt

wMetrics = csv.writer(open("../data/metrics.csv", "w"))
wRuntimes = csv.writer(open("../data/runtimes.csv", "w"))
wRuntimes.writerow(["graphname", "runtime"])

# Get the header, outputting it to the metrics csv file
p = subprocess.Popen(["../bin/MQLib", "-mh"], stdout=subprocess.PIPE)
out, err = p.communicate()
metrics = out.strip().split(",")
wMetrics.writerow(["graphname"] + metrics)

sdb, dom = CloudSetup.setup_sdb_domain("mqlib-metrics")
rs = dom.select('select * from `mqlib-metrics`')
for result in rs:
    timestamp = result["timestamp"].strip()
    graphname = result["graphname"].strip()
    output = result["output"].strip().split(",")
    d = {metrics[i]: float(output[i]) for i in range(len(metrics))}
    rt = runtime(d)
    wMetrics.writerow([graphname] + output)
    wRuntimes.writerow([graphname, rt])
