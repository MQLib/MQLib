import csv
import math
import subprocess
import CloudSetup

wMetrics = csv.writer(open("../data/metrics.csv", "w"))

# Get the header, outputting it to the metrics csv file
p = subprocess.Popen(["../bin/MQLib", "-mh"], stdout=subprocess.PIPE)
out, err = p.communicate()
metrics = out.decode("utf-8").strip().split(",")
wMetrics.writerow(["graphname"] + metrics)

sdb, dom = CloudSetup.setup_sdb_domain("mqlib-metrics")
rs = dom.select('select * from `mqlib-metrics`')
for result in rs:
    timestamp = result["timestamp"].strip()
    graphname = result["graphname"].strip()
    output = result["output"].strip().split(",")
    wMetrics.writerow([graphname] + output)
