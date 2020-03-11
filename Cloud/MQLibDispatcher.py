import boto
import csv
import multiprocessing
import os
import paramiko
import random
import subprocess
import sys
import time

import CloudSetup

def create_instances(tags,verbose=True):
    """
    Simply create an instance for each tag. Uses multiprocessing to create them
    in parallel.
    """

    if verbose: print("Launching instances... ")
    procs = []
    returnInfo = multiprocessing.Queue()
    for tag in tags:
        if verbose: print("  Launching %s ..."%tag)
        proc = multiprocessing.Process(target=CloudSetup.launch_instance,
                                         args=(),
                                       kwargs={"tag":tag,"wait":True,
                                               "returnInfo":returnInfo})
        procs.append(proc)
        proc.start()

    # Wait for all instances to be launched
    for proc in procs:
        proc.join()

    numStarted = 0
    try:
        while True:
            returnInfo.get(False)
            numStarted += 1
    except:
        pass  # Queue is empty
    if numStarted != len(tags):
        print("Exiting because", numStarted, "instances started out of", len(tags))
        exit(0)

    if verbose: 
        print(" All instances launched, but not necessarily ready for")
        print(" SSH though - check AWS console to get a better idea.")
        print(" Hit [RETURN] to proceed to connection attempt.")
        if sys.version_info > (3, 0):
            input()
        else:
            raw_input()

def connect_instances(tags,verbose=True):
    """
    Connect to the instances. Returns, for every tag, an instance handle and
    a cmdshell handle through which we can execute commands and send and
    receive files.
    Returns:
      insts        Dictionary of tag -> instance
      cmds         Dictionary of tag -> cmdshell
    """

    insts = {}
    cmds  = {}
    for tag in tags:
        cmds[tag] = None
    all_done = False
    while not all_done:
        print("Beginning round of connection attempts...")
        all_done = True
        for tag in tags:
            if cmds[tag] == None:
                if verbose: print("  %s"%tag)
                # Test if connected
                try:
                    insts[tag], cmds[tag] = CloudSetup.connect_instance(tag=tag)
                    cmds[tag].run("ls")
                except:
                    all_done = False
                    cmds[tag] = None

    if verbose: 
        print(" All connections established, hit [RETURN]")
        if sys.version_info > (3, 0):
            input()
        else:
            raw_input()
    return insts, cmds


def setup_instances(tags,cmds,verbose=True):
    """
    Install dependencies and build MQLib so it is ready for the run. This takes
    a while so after copying the files we just fire off a script.
    multiprocessing didn't work nicely with this...
    """

    if verbose: print("Copying keys, etc. to all instances, running INSTALL... ")

    # Locate the .boto file
    if os.path.exists(".boto"):
        botoloc = ".boto"
    elif os.path.exists(os.path.expanduser("~/.boto")):
        botoloc = os.path.expanduser("~/.boto")
    else:
        print("Could not locate .boto file")
        exit(1)

    for tag in tags:
        print("    Copying files to %s ..."%(tag))
        f = cmds[tag].open_sftp()
        f.put("INSTALL.sh","INSTALL.sh")  # The install script
        f.put("INSTALL.py", "INSTALL.py")  # Python script that hits INSTALL.sh
        f.put(botoloc, ".boto")  # The AWS key
        f.close()

        print("    Launching INSTALL.py on %s ..."%(tag))
        cmds[tag].run("chmod +x INSTALL.sh")  # Make script executable        
        stdin, stdout, stderr = cmds[tag]._ssh_client.exec_command("python INSTALL.py " + sys.argv[3])  # Not blocking

    print("    Waiting for INSTALL.py to complete on all machines")
    while True:
        time.sleep(10)
        done = [cmds[tag].run("ls MQLib/bin/")[1].decode("utf-8").find("MQLib") >= 0 for tag in tags]
        print("      - Installation complete on", sum(done), "/", len(done), "boxes")
        if sum(done) == len(done):
            break

    if verbose:
        print(" Hit [RETURN] when ready to proceed.")
        if sys.version_info > (3, 0):
            input()
        else:
            raw_input()

def split_graphs(tags):
    """
    Split up graphs between the machines. Deterministic split so won't change
    as long as the graphs uploaded to S3 don't change.
    """
    mach_graphs = {}
    for tag in tags:
        mach_graphs[tag] = []
    # Load in the graphs -- load from file for a validation run and otherwise
    # load from s3.
    if sys.argv[1] == "FULL":
        try:
            instance_file = "../data/instances.txt"
            all_graphs = [x.strip() for x in open(instance_file, "r")]
            # Reproducibly shuffle the graphs to distribute the workload
            random.seed(144)
            random.shuffle(all_graphs)
            print("Loaded", len(all_graphs), "instances from", instance_file)
        except Exception as e:
            print("You need to create a", instance_file, "file with one")
            print("  graph per line for a full run.")
            exit(1)
    elif sys.argv[1] == "METRICS" and os.path.exists("metric_files.txt"):
        all_graphs = sorted([x.strip() for x in open("metric_files.txt", "r")])
        print("Loaded", len(all_graphs), "instances from metric_files.txt")
    else:
        conn = boto.connect_s3(anon=True)
        b = conn.get_bucket("mqlibinstances", validate=False)
        all_graphs = sorted([key.name for key in b.list()])

    # Split them up (add graphs to each node first in tag order and then in
    # reverse tag order to balance the workload across nodes for the full run)
    num_graphs = len(all_graphs)
    cur_graph_ind = 0
    while cur_graph_ind < num_graphs:
        for tag in tags:
            mach_graphs[tag].append(all_graphs[cur_graph_ind])
            cur_graph_ind += 1
            if cur_graph_ind == num_graphs:
                for tag in tags:
                    mach_graphs[tag] = mach_graphs[tag][::-1]  # Increasing
                return mach_graphs
        for tag in reversed(tags):
            mach_graphs[tag].append(all_graphs[cur_graph_ind])
            cur_graph_ind += 1
            if cur_graph_ind == num_graphs:
                for tag in tags:
                    mach_graphs[tag] = mach_graphs[tag][::-1]  # Increasing
                return mach_graphs

def dispatch_and_run(tags,cmds,mach_graphs,verbose=True):
    """
    Kick things off by copying a file listing all graphs this instance
    must run and run the MQLibMaster
    """
    # Write out and copy to instances
    if verbose: print("Copying graphs and heuristics files and starting run... ")
    for tag in tags:
        if verbose: print(" %s"%tag)

        with open("/tmp/GRAPH_FILE_%s"%tag,"w") as fp:
            for g in mach_graphs[tag]:
                fp.write(g + "\n")

        f = cmds[tag].open_sftp()
        f.put("/tmp/GRAPH_FILE_%s"%tag, "MQLib/Cloud/GRAPH_FILE")
        if sys.argv[1] == "FULL":
            f.put("../data/heuristics.txt", "MQLib/Cloud/HEUR_FILE")
        f.close()

        if sys.argv[1] == "FULL":
            cmds[tag]._ssh_client.exec_command("cd MQLib/Cloud; nohup python MQLibMaster.py FULL " + tag + " " + sys.argv[4] + " " + sys.argv[5] + " " + sys.argv[6] + " " + sys.argv[7] + " &> screen_output.txt &")
        else:
            cmds[tag]._ssh_client.exec_command("cd MQLib/Cloud; nohup python MQLibMaster.py " + sys.argv[1] + " " + tag + " &> screen_output.txt &")
    if verbose: print("\n  Computation started on all machines")

def run_dispatch():
    """
    Setup machines, run jobs, monitor, then tear them down again.
    """

    if len(sys.argv) < 4 or not sys.argv[1] in ["METRICS", "FULL"] or not sys.argv[2].isdigit() or sys.argv[3].find(".git") < 0 or (sys.argv[1] == "FULL" and (len(sys.argv) < 8 or not sys.argv[4].isdigit() or not all([x.isdigit() for x in sys.argv[5].split("_")]) or not sys.argv[6].lstrip("-").isdigit() or not sys.argv[7].lstrip("-").isdigit())):
        print("Usage:\n  python MQLibDispatcher.py METRICS #NODE http://link/to/git/repo.git [nocreate] [nodispatch] [verbose]\n    [[or]]\n  python MQLibDispatcher.py FULL #NODE http://link/to/git/repo.git #ITERFORBASELINE SEEDS_SEPARATED_BY_UNDERSCORES MINSECONDS MAXSECONDS [nocreate] [nodispatch] [verbose]")
        exit(1)
    NUM_MACHINES     = int(sys.argv[2])
    tags = ["mqlibtest%d"%i for i in range(NUM_MACHINES)]
    CREATE_MACH      = not ("nocreate" in sys.argv)
    DISPATCH_AND_RUN = not ("nodispatch" in sys.argv)
    VERBOSE          = "verbose" in sys.argv

    # Validate that required files exist
    if not os.path.exists(".boto") and not os.path.exists(os.path.expanduser("~/.boto")):
        print("Please create a .boto file containing your AWS credentials as")
        print(" described in file Cloud/README.md. Store this file either in")
        print(" the Cloud folder or in your home directory.")
        exit(1)
    if not os.path.exists("INSTALL.sh") or not os.path.exists("INSTALL.py"):
        print("Please run this script from the Cloud directory.")
        exit(1)
    if sys.argv[1] == "FULL":
        instance_file = "../data/instances.txt"
        if not os.path.exists(instance_file):
            print("You need to create a", instance_file, "file with one")
            print("  instance per line for a full run.")
            exit(1)
        heuristic_file = "../data/heuristics.txt"
        if not os.path.exists(heuristic_file):
            print("You need to create a", heuristic_file, "file with one")
            print("  heuristic to be tested per line.")
            exit(1)

    # Split up graphs between the machines.
    mach_graphs = split_graphs(tags)
    
    # Setup security group and key pair (these are no-ops if done before)
    CloudSetup.create_security_group()
    CloudSetup.create_keypair()

    # Do some work upfront (clean out ~/.ssh/known_hosts and wait until all
    # shutting down nodes are shut down).
    CloudSetup.clean_known_hosts()
    CloudSetup.wait_for_shutdown()

    # Create instances if desired
    if CREATE_MACH: create_instances(tags, VERBOSE)

    # Connect to all the instances
    insts, cmds = connect_instances(tags,VERBOSE)

    # Set them up (if desired)
    if CREATE_MACH: setup_instances(tags,cmds,VERBOSE)
        
    # Send out jobs and start machines working (if desired)
    if DISPATCH_AND_RUN: dispatch_and_run(tags,cmds,mach_graphs,VERBOSE)

    print("")
    print("All dispatcher tasks successfully completed.")

if __name__ == "__main__":
    run_dispatch()
