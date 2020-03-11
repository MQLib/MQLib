import boto
import datetime
import os
import shutil
import subprocess
import sys
import zipfile
import CloudSetup

def get_graph(graph):
    """
    Abstracts out the process of pulling the graph from main test loop
    """
    print("Acquiring graph %s..."%graph)
    conn = boto.connect_s3(anon=True)
    b = conn.get_bucket("mqlibinstances", validate=False)
    k = boto.s3.key.Key(b)
    k.key = graph
    k.get_contents_to_filename("curgraph.zip")
    with zipfile.ZipFile("curgraph.zip", "r") as z:
        for name in z.namelist():
            outloc = z.extract(name)
            shutil.move(outloc, "curgraph")
    print("    done!")
    sys.stdout.flush()

def chunks(l, n):
    """
    Break up a string l into chunks of length at most n, returning a list of
    all the chunks.
    """

    if n < 1:
        n = 1
    return [l[i:i + n] for i in range(0, len(l), n)]

def filterOutput(output, consec, ending):
    """
    Filter the output of a heuristic by only reporting every "consec" new best solution,
    as well as the last "ending" solutions.
    """
    intro = output.split("[")[0]
    toFilter = output.strip().split("[")[1][:-1].split(";")
    filterIdx = range(0, len(toFilter)-ending, consec) + range(len(toFilter)-ending, len(toFilter))
    filtered = [toFilter[x] for x in filterIdx]
    return intro + "[" + ";".join(filtered) + "]"

def run_tasks():
    """
    Run all the tasks allocated to this worker.

    Graphs are assigned to workers, then that worker runs all heuristics
    multiple times on that graph. The graph is then pulled from the S3
    datastore and the results are inserted into a SimpleDB domain, with
        key: [graphname]-[heuristic]-[run]
        value: {"graphname" : graph_name,
                "heuristic" : heuristic,  # Only for run type "FULL"
                "run"       : run_id,  # Only for run type "FULL"
                "timestamp" : timestamp,
                "ouput"     : output_from_test  # Can be chunked (see below)}
    When deciding what to run, a check will be made for existing runs
    in SimpleDB to avoid duplicating work, making restarts easier
    """

    # Validate command-line arguments
    if not (sys.argv[1] == "METRICS" and len(sys.argv) == 3) and not (sys.argv[1] == "FULL" and len(sys.argv) == 7 and sys.argv[3].isdigit() and all([x.isdigit() for x in sys.argv[4].split("_")]) and sys.argv[5].lstrip("-").isdigit() and sys.argv[6].lstrip("-").isdigit()):
        print("Usage:\n  python MQLibRunner.py METRICS tag\n    [[or]]\n  python MQLibRunner.py FULL tag #ITERFORBASELINE SEEDS_SEPARATED_BY_UNDERSCORES MINSECONDS MAXSECONDS")
        exit(1)
    if sys.argv[1] == "METRICS":
        SDB_DOMAIN = "mqlib-metrics"
    else:
        SDB_DOMAIN = "mqlib-domain2"

    # Graphs I need to run will be in a local file called GRAPH_FILE
    with open("GRAPH_FILE","r") as fp:
        lines = fp.readlines()
        if sys.argv[1] == "FULL":
            GRAPHS = [line.strip() for line in lines]
        else:
            GRAPHS = [line.strip() for line in lines]
    print(GRAPHS)

    # Heuristics I need to run - shared by all, as all workers run all
    # heuristics (possibly a good thing for comparability too)
    if sys.argv[1] == "FULL":
        with open("HEUR_FILE","r") as fp:
            HEURS = [line.strip() for line in fp.readlines()]
        print(HEURS)
        SEEDS = [int(x) for x in sys.argv[4].strip().split("_")]
        print(SEEDS)
        sys.stdout.flush()

    # Get ready for SimpleDB comms
    sdb, dom = CloudSetup.setup_sdb_domain(SDB_DOMAIN)
    
    # Log progress to a file so we can be monitored easily
    log_fp = open("PROGRESS","w")

    for gidx in range(len(GRAPHS)):
        graph = GRAPHS[gidx]
        if sys.argv[1] == "FULL":
            # First, acquire or build the graph and place in local folder
            # as 'curgraph'
            get_graph(graph)
            log_fp.write(graph + ", " + str(datetime.datetime.now()))
            log_fp.flush()

            # Run the baseline to figure out the runtime limit
            torun = ["timeout", "10000", "../bin/MQLib", "-fM",
                     "curgraph", "-h", "BASELINE", "-r",
                     sys.argv[3], "-s", "144"]
            print(torun)
            sys.stdout.flush()
            p = subprocess.Popen(torun, stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT)
            baseline_output = p.stdout.read()
            runtime = baseline_output.split(",")[4]

            if int(sys.argv[5]) >= 0:
                runtime = str(max(int(sys.argv[5]), float(runtime)))
            if int(sys.argv[6]) >= 0:
                runtime = str(min(int(sys.argv[6]), float(runtime)))

            print("******** Baseline run:***")
            print(baseline_output)
            print("")
            print("******** Extracted runtime limit:")
            print(runtime)
            
            for heur in HEURS:
                # See how much work we've done for this heuristic before
                reps_to_do = {}
                for seed in SEEDS:
                    reps_to_do[seed] = True
                query = 'SELECT * FROM `%s` WHERE heuristic="%s" AND graphname="%s"'%(SDB_DOMAIN,heur,graph)
                #print(query)
                query_results = dom.select(query)
                for q_r in query_results:
                    #print(q_r)
                    reps_to_do[int(q_r["run"])] = False

                # Do the reps we haven't done yet
                for seed in SEEDS:
                    if not reps_to_do[seed]:
                        print("Skipping rep %d for %s - %s"%(seed,graph,heur))
                        continue

                    # Run the program
                    print("Running rep %d for %s - %s"%(seed,graph,heur))
                    torun = ["timeout", "10000", "../bin/MQLib", "-fM",
                             "curgraph", "-h", heur, "-r", runtime, "-s",
                             str(seed)]
                    print(torun)
                    sys.stdout.flush()
                    p = subprocess.Popen(torun, stdout=subprocess.PIPE,
                                         stderr=subprocess.STDOUT)
                    mqlib_output = p.stdout.read()

                    key = "%s-%s-%d"%(graph,heur,seed)
                    value = {   "graphname" : graph,
                                "heuristic" : heur,
                                "run"       : str(seed),
                                "timestamp" : str(datetime.datetime.now())
                            }

                    # Because mqlib_output might exceed the limit of length of
                    # 1024 for sdb attributes, we'll split it up into
                    # 1000-character chunks and add each chunk as attributes
                    # output000, output001, ...
                    parts = chunks(mqlib_output, 1000)

                    # If the output is too large, filter the solutions to take
                    # every "consec" reported value, as well as the last 10.
                    # Keep increasing "consec" until the output is of acceptable
                    # size (aka won't cause an error when loaded into the SimpleDB).
                    consec = 1
                    while len(parts) >= 55:
                        consec += 1
                        parts = chunks(filterOutput(mqlib_output, consec, 10), 1000)

                    for num in range(len(parts)):
                        attr = "output" + str(num).zfill(3)
                        value[attr] = parts[num]

                    dom.put_attributes(key, value)
        else:
            # Metrics run
            query = 'SELECT * FROM `%s` WHERE graphname="%s"'%(SDB_DOMAIN,graph)
            query_results = dom.select(query)
            if next(query_results, None):
                print("Skipping %s" % graph)
                continue

            # First, acquire or build the graph and place in local folder
            # as 'curgraph'
            get_graph(graph)
            log_fp.write(graph + ", " + str(datetime.datetime.now()))
            log_fp.flush()

            # No need to set seed on metrics run because the metrics code always
            # sets the seed to 0 before running.
            torun = ["../bin/MQLib", "-fM", "curgraph", "-m"]
            sys.stdout.flush()
            p = subprocess.Popen(torun, stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT)
            mqlib_output = p.stdout.read()
            value = {"graphname" : graph,
                     "output"    : mqlib_output,
                     "timestamp" : str(datetime.datetime.now())
                     }
            dom.put_attributes(graph, value)

        log_fp.write(", " + str(datetime.datetime.now()) + "\n")
        log_fp.flush()

    log_fp.close()

    # Self-terminate at completion
    CloudSetup.terminate_instance(sys.argv[2])

if __name__ == "__main__":
    run_tasks()
