# These are the run-once tasks to get AWS ready-for-use, in addition to utilities
# used by our scripts to access AWS.
import csv
import os
import time
import traceback
import boto.ec2
import boto.manage.cmdshell
import boto.sdb

SSH_FOLDER = os.path.expanduser("~/.ssh/")
INST_TYPE  = "m3.medium"
AWS_REGION = "us-east-1"  # US East (Virginia)
AMI_NAME   = "ami-08531832555a14fe2"    # US East (Virginia) Ubuntu 14.04 LTS (HVM)

def create_keypair(key_name="mqlibtest"):
    """
    Create the public-key crypto pair so we can log in to our new instances.
    AWS stores the public key under a name we provide, we need to save the
    private key ourselves.
    """
    if os.path.isfile(SSH_FOLDER + key_name + ".pem"):
        return  # Key already created
    ec2 = boto.ec2.connect_to_region(AWS_REGION)
    key = ec2.create_key_pair(key_name)
    key.material = key.material.encode()  # https://github.com/boto/boto/issues/3782
    key.save(SSH_FOLDER)

def create_security_group(group_name="mqlibtest"):
    """
    Instances are pretty locked down by default. We can assign them to
    security groups to give access rights. This group allows us to access
    our instances over SSH
    """
    ec2 = boto.ec2.connect_to_region(AWS_REGION)
    for g in ec2.get_all_security_groups():
        if g.name == group_name:
            return  # We already have this group setup
    group = ec2.create_security_group(
                group_name, "MQLib SSH access group")
    group.authorize("tcp", 22, 22, "0.0.0.0/0")  # SSH is on port 22, all IPs
    print("Created new security group")

def launch_instance(tag="mqlibtest1",key_name="mqlibtest",group_name="mqlibtest",
                    wait=True, returnInfo=None):
    """
    Launch a testing instance. Doesn't actually attempt to connect as
    it can take quite a while between 'running' and connectability
    """
    ec2 = boto.ec2.connect_to_region(AWS_REGION)
    failures = 0
    max_failures = 10
    while True:
        try:
            reservation = ec2.run_instances(AMI_NAME,
                                            key_name=key_name,
                                            security_groups=[group_name],
                                            instance_type=INST_TYPE,
                                            user_data=None)
            break
        except Exception as err:
            # Failed to get instance; wait 15 seconds and then try again (up to
            # 10 total times)
            errorText = str(err)
            if errorText.find("Not authorized for images") >= 0:
                print("**************************************")
                print("* Error from AWS suggests that the AMI code in")
                print("* CloudSetup.py is deprecated. Please go to")
                print("* https://aws.amazon.com/marketplace/ and search for")
                print("* \"Ubuntu server lts hvm\", selecting the most recent")
                print("* version. Click \"Continue to Subscribe\",")
                print("* \"Accept Terms\", \"Continue to Configuration\"")
                print("* and then copy the AMI ID for the US East region.")
                print("* Copy that to the AMI_NAME value in CloudSetup.py")
                print("* and re-run.")
                print("***************************************")
                print("* (Full text of error):")
                print(errorText)
                print("***************************************")
                return None
            elif errorText.find("accept terms and subscribe") >= 0:
                print("**************************************")
                print("* Error from AWS suggests that you have never used this")
                print("* AMI before and need to accept its terms and")
                print("* subscribe to it. Please follow the link in the below")
                print("* error text. Click \"Continue\", \"Manual Launch\",")
                print("* and \"Accept Terms\". After receiving email")
                print("* confirmation, you can re-run the code.")
                print("**************************************")
                print("* (Full text of error):")
                print(errorText)
                print("**************************************")
                return None
            failures += 1
            if failures == max_failures:
                print("**************************************")
                print("* Maximum number of instance launch failures reached.")
                print("* (Full text of error):")
                print(errorText)
                print("**************************************")
                return None
            print("    ** ec2.run_instances failed for tag", tag, "; waiting 15")
            print("    ** seconds and then trying again...")
            time.sleep(15)

    time.sleep(5)  # Slow things down -- they're never running super fast anyway
    instance = reservation.instances[0]
    time.sleep(5)  # Slow things down -- they're never running super fast anyway
    instance.add_tag("mqlibname",tag)
    time.sleep(5)  # Slow things down -- they're never running super fast anyway

    if wait:
        print("    Instance requested, waiting for 'running' for tag", tag)
        while instance.state != "running":
            print("    %s ..."%tag)
            time.sleep(5)
            try:
                instance.update()
            except EC2ResponseError as e:
                print("******************")
                print("Error caught in instance.update():")
                print(e.strerror)
                print("******************")
        print("    %s done!"%tag)
    if returnInfo:
        returnInfo.put(tag)
    return instance

def get_instance(tag="mqlibtest1"):
    """
    Get instance by tag
    """
    ec2 = boto.ec2.connect_to_region(AWS_REGION)
    reservations = ec2.get_all_instances()
    for res in reservations:
        for inst in res.instances:
            if "mqlibname" in inst.tags.keys():
                if inst.tags["mqlibname"] == tag and inst.state == "running":
                    #print("Found %s"%tag)
                    return inst
    print("Couldn't find instance")
    return None

def connect_instance(tag="mqlibtest1",key_name="mqlibtest"):
    """
    Connect to a running instance using a tag
    """
    inst = get_instance(tag)
    cmd = boto.manage.cmdshell.sshclient_from_instance(inst,
                SSH_FOLDER+key_name+".pem",
                user_name="ubuntu")
    return inst, cmd

def terminate_instance(tag="mqlibtest1"):
    inst = get_instance(tag)
    inst.terminate()

###############################################################################

def setup_sdb_domain(domain_name="mqlib-domain"):
    sdb = boto.sdb.connect_to_region(AWS_REGION)
    # Only create if it doesn't exist already
    try:
        dom = sdb.get_domain(domain_name,validate=True)
    except:
        # Doesn't exist yet
        dom = sdb.create_domain(domain_name)
    return sdb, dom

def delete_sdb_domain(domain_name="mqlib-domain"):
    sdb, dom = setup_sdb_domain(domain_name)
    sdb.delete_domain(domain_name)

def dump_sdb_domain(domain_name="mqlib-domain"):
    sdb, dom = setup_sdb_domain(domain_name)
    rs = dom.select('select * from `' + domain_name + '`')
    for j in rs:
        print(j)

def get_sdb_domain_size(domain_name="mqlib-domain"):
    sdb, dom = setup_sdb_domain(domain_name)
    rs = dom.select('select count(*) from `' + domain_name + '`')
    ct = 0
    for res in rs:
        ct += int(res[u'Count'])
    print("Size of", domain_name, ":", ct)

def copy_full_run(from_domain, to_domain):
    sdbFrom, domFrom = setup_sdb_domain(from_domain)
    sdbTo, domTo = setup_sdb_domain(to_domain)
    rs = domFrom.select('select * from `' + from_domain + '`')
    
    numPut = 0
    values = {}
    for x in rs:
        key = x["graphname"] + "-" + x["heuristic"] + "-" + x["run"]
        values[key] = x
        if len(values) == 25:
            domTo.batch_put_attributes(values)
            numPut += len(values.keys())
            print(numPut)
            values = {}
    if len(values.keys()) > 0:
        domTo.batch_put_attributes(values)
    print("Done")

def remove_graphs_with_names(domain_name, graphs):
    sdb, dom = setup_sdb_domain(domain_name)
    values = {}
    for idx in range(len(graphs)):
        graph = graphs[idx]
        print("Removing runs from", graph, "(", idx+1, "/", len(graphs), ")")
        rs = dom.select('select * from `' + domain_name + '` where graphname="' + graph + '"')
        for x in rs:
            key = x["graphname"] + "-" + x["heuristic"] + "-" + x["run"]
            if key in values:
                print("OVERLAP!!")
            values[key] = None
            if len(values.keys()) == 25:
                dom.batch_delete_attributes(values)
                values = {}
        if len(values.keys()) > 0:
            dom.batch_delete_attributes(values)

def remove_graphs_with_file(domain_name, filename):
    reader = csv.reader(open(filename, "rU"))
    header = reader.next()
    if header != ["graphname", "runtime"]:
        print("Illegal header in remove_graphs_with_file")
        exit(0)
    graphs = [x[0] for x in reader]
    remove_graphs_with_names(domain_name, graphs)

# Several cleanup tasks to make starting a cluster less annoying:
def clean_known_hosts():
    with open(SSH_FOLDER + "known_hosts", "rU") as fp:
        lines = fp.readlines()
        filtered = [x for x in lines if x.find("ec2-") != 0]
    with open(SSH_FOLDER + "known_hosts", "w") as fp:
        for line in filtered:
            fp.write(line)
    print("Removed", len(lines) - len(filtered), "lines from ~/.ssh/known_hosts")

def get_num_running():
    ec2 = boto.ec2.connect_to_region(AWS_REGION)
    reservations = ec2.get_all_instances()
    num_shutting_down = 0
    num_pending_running = 0
    num_stop = 0
    num_terminate = 0
    for res in reservations:
        for inst in res.instances:
            if inst.state == "shutting-down":
                num_shutting_down += 1
            elif inst.state in ["pending", "running"]:
                num_pending_running += 1
            elif inst.state in ["stopping", "stopped"]:
                num_stop += 1
            elif inst.state == "terminated":
                num_terminate += 1
    return (num_shutting_down, num_pending_running, num_stop, num_terminate)

def print_num_running():
    nr = get_num_running()
    print("Number Shutting Down:", nr[0])
    print("Number Pending or Running:", nr[1])
    print("Number Stopping or Stopped:", nr[2])
    print("Number Terminated:", nr[3])

# Wait for all the EC2 nodes that are in shutting-down to go to status 
def wait_for_shutdown():
    while True:
        n_shut_down, n_pend_run, n_stop, n_terminate = get_num_running()
        if n_shut_down == 0:
            print("No nodes shutting down")
            return
        else:
            print(n_shut_down, "instance(s) still shutting down and", n_pend_run, "pending/running; waiting")
            time.sleep(5.0)
