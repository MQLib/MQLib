# Reproducible Parallel Computation with Amazon Web Services

This guide describes how to perform reproducible experimentation using Amazon Web Services (AWS).

## Setting Up Required Software

You will need to install several pieces of software on your system to interface with AWS:

 * Python version 2.7 should be installed
 * If you do not already have the Python package manager `pip` installed, install it by following the instructions [on the pip website](https://pip.pypa.io/en/latest/installing.html).
 * Install required python packages with `sudo pip install paramiko` and `sudo pip install boto`. The `paramiko` package dependencies `ecdsa` and `pycrypto` should be automatically installed; if they are not follow the additional [paramiko installation instructions](http://www.paramiko.org/installing.html).

## Configuring an Amazon Web Services Account

The first step is to set up an Amazon Web Services account:

 * Register for an account at [aws.amazon.com](https://aws.amazon.com) (registration will require credit card information).
 * Log in to the AWS console at [console.aws.amazon.com](https://console.aws.amazon.com).
 * Navigate to "Identity & Access Management"; click "Users" and "Create New Users", creating one new user with "Generate an access key for each user" checked.
 * Click "Show User Security Credentials", which should show a "Access Key ID" value and a "Secret Access Key" value. Use these values to create a file located at `~/.boto` with the following format (filling in the *'s with the value from the user security credentials).
```
[Credentials]
aws_access_key_id = **************
aws_secret_access_key = **************
```
 * From the console, navigate to "Identity & Access Management"; click "Users" and select the new user. Under "Inline Policies", click the text "click here" to add a new inline policy (you may need to click the arrow in the "Inline Policies" heading to display this text). Click "Custom Policy" and "Select", entering any policy name desired and the following policy (which grants the user full access to the account):
```
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": "*",
      "Resource": "*"
    }
  ]
}
```

## Generating Problem Instance Metrics

The first type of cloud-based computation is computing the metrics for each problem instance stored in the publicly accessible `mqlibinstances` S3 bucket, creating a file [data/metrics.csv](../data/metrics.csv) containing the metrics for each problem instance and a file [data/runtimes.csv](../data/runtimes.csv) containing the runtime limits for each problem instance. Versions of these files containing the metrics and runtime limit formula from [Dunning et al. (2015)] are checked into the repository, so there is no need to re-compute these values unless changing the metrics or termination criteria.

### Computing metric values on EC2

To compute problem instance metrics, the first steps are to implement the new metric and to check it into a forked version of the MQLib repository, as detailed in the [MQLib guide](../README.md) under the "Contributing a New Problem Instance Metric" heading.

The [Cloud/MQLibDispatcher.py](MQLibDispatcher.py) script enables you to perform a distributed computation of the metric values using Amazon EC2. The simplest way to perform this computation is to run this script from the `Cloud` directory using the command line, specifying the `METRICS` task, the number of EC2 nodes to launch, and the git repository to use on each of the EC2 nodes. In general this will be a forked version of the repository located at `https://github.com/USERNAME/MQLib.git`, where `USERNAME` is your GitHub username; refer to the [MQLib guide](../README.md) for detailed instructions on how to determine the clone URL for the forked repository. If you instead wanted to use the main repository, you could specify `https://github.com/MQLib/MQLib.git`.

```
python MQLibDispatcher.py METRICS 10 https://github.com/USERNAME/MQLib.git
```

This will produce output as the dispatcher launches new nodes, connects to the nodes, sets up the nodes for the run, distributes the work for the run (graphs to have their metrics computed in this case), and starts the run on each node. Once the launch script exits, all nodes will be running, and they will terminate automatically once they have completed all their assigned work. See the "Monitoring Cloud Runs" section below for details on how to monitor a run that is in progress.

By default, this run will compute the metrics for all the problem instances stored in the publicly accessible S3 bucket `mqlibinstances`. To specify a subset of the instances in this bucket, create a text file called `metric_files.txt` in the `Cloud` folder with the name of each specified instance, one per line. To see the list of all available instances, you can run `python scripts/list_graphs.py`. An example file would be:

```
p3000_2.zip
imgseg_43070.zip
imgseg_43083.zip
```

By default an account is limited to 20 simultaneous EC2 nodes of a given type, though an increase can be requested on the AWS Console by selecting the "US East" region, selecting EC2, and then selecting "Limits" on the navigation pane on the left, requesting a limit increase for the "m3.medium" instance type.

### Downloading metric values and computing termination criteria

Once the run is completed, its full output is stored in a SimpleDB (NoSQL database) on Amazon Web Services called `mqlib-metrics`. To download these metrics to `data/metrics.csv` and to compute the runtime limits from these metrics, storing the results in `data/runtimes.csv`, you can run the following from the `Cloud` directory:

```
python grabGraphInfo.py
```

The script computes the runtime limit using the function from [Dunning et al. (2015)], but the function `runtime` in the script could be replaced with any other function of the metrics, enabling other runtime limits to be used.

## Testing Heuristics

The main type of cloud-based computation is running a specified set of heuristics over a specified set of problem instances, using instance-specific runtime limits as the termination criteria. This process will generate a raw datafile containing the elapsed runtime and solution value of every new best solution identified during the solution process in `data/results.csv` (or a similar name to avoid overwriting an existing file). A zipped version of this dataset containing the results for all combinations of the 37 heuristics and 3,398 problem instances considered in [Dunning et al. (2015)] is provided at [data/results.zip](../data/results.zip).

### Testing heuristics on EC2

To test a specified set of heuristics in a distributed manner on EC2, you can again run the [Cloud/MQLibDispatcher.py](MQLibDispatcher.py) script from the `Cloud` directory using the command line, specifying the number of nodes to launch on EC2 and the git repository to use on each of the EC2 nodes. Again, in general this will be a forked version of the repository located at `https://github.com/USERNAME/MQLib.git`, where `USERNAME` is your GitHub username; refer to the [MQLib guide](../README.md) for detailed instructions on how to determine the clone URL for the forked repository. If you instead wanted to use the main repository, you could specify `https://github.com/MQLib/MQLib.git`.

```
python MQLibDispatcher.py FULL 10 https://github.com/USERNAME/MQLib.git
```

This script requires the following files to be setup to function properly:

1. [data/runtimes.csv](../data/runtimes.csv): This file specifies the problem instances to test and the runtime limit for each. The version checked into the repository contains all the files in the public `mqlibinstances` S3 bucket with the runtime limits from [Dunning et al. (2015)], and the "Downloading metric values and computing termination criteria" section above in this document describes how to compute new runtime limits using the [Cloud/grabGraphInfo.py](grabGraphInfo.py) script.
2. [data/heuristics.txt](../data/heuristics.txt): This file specifies the heuristic codes to be tested, with one heuristic per line. As stated in the [Using the MQLib guide](../bin/README.md), you can view all available heuristic codes by running `./bin/MQLib -l`, and the list of all 37 heuristics tested in [Dunning et al. (2015)] are checked into the repository. To test a new heuristic, this file should contain a single line with that heuristic's code.

To prevent infinite looping scripts from derailing the automatic testing, each heuristic is given a 1-hour time limit to perform testing on an instance (this includes the time required to read in the instance).

### Downloading heuristic results

Once the run is completed, its full output (the set of all new best solutions for all heuristic/graph pairs) is stored in a SimpleDB (NoSQL database) on Amazon Web Services called `mqlib-domain`. To download these results to `data/results.csv` and a list of errors encountered in the run to `data/errors.txt` (or numbered versions of these files if they already exist), you can run the following from the `Cloud` directory:

```
python grabFullRuns.py ../data/runtimes.csv
```

The argument to this script is the csv file of instance-specific runtime limits that was used for the run.

## Monitoring and Debugging Cloud Runs

There are two major places where a cloud run can fail: during node setup and during the runs themselves.

### Debugging Node Setup Failures

The first time you launch a run on the cloud, it should fail with a message saying that you need to accept the terms and conditions of the Amazon Machine Image (AMI) that we use on EC2. In this case, simply follow the outputted instructions and re-run.

Periodically Amazon deprecates AMIs and releases newer versions of them, causing node launches to fail. In this case, please follow the outputted instructions to find the updated AMI and to update the `AMI_NAME` variable in [Cloud/CloudSetup.py](CloudSetup.py). When re-running you should be prompted to accept the terms and conditions of the new AMI, and after that you should be able to successfully run. Please consider creating a GitHub pull request with the updated `AMI_NAME` variable so others can avoid this error.

Otherwise, the most likely error to be encountered is one in which the `MQLibDispatcher.py` script outputs the following ad infinitum (where `xx` is the number of machines being run on; there is probably trouble if this message appears at least 100 times):

```
      - Installation complete on 0 / xx boxes
      - Installation complete on 0 / xx boxes
      - Installation complete on 0 / xx boxes
      - Installation complete on 0 / xx boxes
      - Installation complete on 0 / xx boxes
      ...
```

This likely indicates an issue in [Cloud/INSTALL.sh](INSTALL.sh), the script that performs the main setup tasks on each EC2 node. To resolve such an issue, either contact the MQLib package maintainers or debug by logging onto a node and looking at the log output generated by [Cloud/INSTALL.sh](INSTALL.sh). The first step to log into a node is to identify that node's domain name. You can do this by logging into the AWS console, selecting "EC2", selecting "Running Instances", selecting an instance, and reading the value under "Public DNS" (we will call this `DNS` in the command that follows). Then you can log into the instance on the command line with:

```
ssh -i ~/.ssh/mqlibtest.pem ubuntu@DNS
```

Here, `~/.ssh/mqlibtest.pem` is a key that was generated when running [Cloud/MQLibDispatcher.py](MQLibDispatcher.py). It is created by the `CloudSetup.create_keypair` function to enable access to the EC2 nodes without a password.

Once you have logged onto a node, there will be files containing output from the setup tasks with the following names:

* `progress_A_1.txt`: Output from updating packages with `sudo apt-get update --yes`.
* `progress_B_1.txt`: Output from installing packages `git`, `build-essentials`, `python-pip`, and `python-paramiko` with `apt-get`.
* `progress_C_1.txt`: Output from installing the `boto` package with `pip`.
* `progress_D_1.txt`: A file indicating that the `rm -rf MQLib` call was successfully completed.
* `progress_E_1.txt`: Output from cloning the indicated git repository.
* `progress_F_1.txt`: Output from building the MQLib.

If these steps failed to produce an `MQLib` executable, then there may be additional files with output from additional setup attempts (numbered 2, 3, ...). From reading this output, you may be able to identify and correct a setup issue, with the most likely culprit being that the GitHub repository was incorrectly specified on the command line. Once you determine the cause of the error, you can terminate all running instances from the AWS console by navigating to "EC2" and "Running Instances", selecting all the instances, right clicking, and selecting "Instance State -> Terminate".

If the issue you were investigating was due to a bug, please consider submitting a GitHub pull request correcting the bug. If you cannot determine the source of the issue, please contact the MQLib team.

### Monitoring and Debugging After Setup

After the [Cloud/MQLibDispatcher.py](MQLibDispatcher.py) script has finished executing, there are two main options for monitoring your processes:

* At the crudest level, you can check the number of running processes, as EC2 nodes initialized by [Cloud/MQLibDispatcher.py](MQLibDispatcher.py) continue running until they have finished all their assigned work, at which point they terminate. One way to check the status of the EC2 nodes is logging onto the AWS console, selecting "EC2" and "Running Instances". Another way is running python interactively from the `Cloud` folder and then running `import CloudSetup` and then `CloudSetup.print_num_running()`.
* To make sure the EC2 instances are making progress, you can check the number of tasks they have completed. Runs of [Cloud/MQLibDispatcher.py](MQLibDispatcher.py) in `METRICS` mode output to the `mqlib-metrics` SimpleDB, while runs in `FULL` mode output to the `mqlib-domain` SimpleDB. You can check the sizes of these databases by running python interactively from the `Cloud` folder, running `import CloudSetup`, and then either running `CloudSetup.get_sdb_domain_size('mqlib-metrics')` or running `CloudSetup.get_sdb_domain_size('mqlib-domain')`.

If there is no progress, then you can debug the run by logging onto one of the EC2 nodes as described in the previous section of this document. Debug output will be available in files `MQLib/Cloud/screen_output.txt` and `MQLib/Cloud/PROGRESS`. Once you are finished debugging, you can terminate all running instances from the AWS console by navigating to "EC2" and "Running Instances", selecting all the instances, right clicking, and selecting "Instance State -> Terminate".

## Internals

Though we do not expect that most researchers will need to modify the internals of the AWS testing infrastructure, we provide details for completeness. As indicated in the previous sections, cloud runs are started with [Cloud/MQLibDispatcher.py](MQLibDispatcher.py). This script performs a number of steps:

1. Validates that the command-line options are properly specified, that the `.boto` credentials file is present, and that the script is run from the `Cloud` directory. In the case of a `FULL` run, verifies that the required files ([data/runtimes.csv](../data/runtimes.csv) and [data/heuristics.txt](../data/heuristics.txt)) are present.
2. Performs several startup tasks. The `CloudSetup.create_security_group` function is used to create an AWS security group (if one has not already been created) that allows SSH access to the EC2 nodes from any IP address on port 22. The `CloudSetup.create_keypair` function is used to create a key for the account's user (if one has not already been created), which is stored in `~/.ssh/mqlibtest.pem`. The `CloudSetup.clean_known_hosts` function is used to remove any EC2 hosts from the `~/.ssh/known_hosts` file, which prevents SSH errors due to hostname collisions in consecutive large runs. Finally, the `CloudSetup.wait_for_shutdown` function waits for all nodes that are currently shutting down to be terminated.
3. Creates the EC2 instances using the `create_instances` function (this step is skipped if the `nocreate` command-line argument is provided). This is a thin wrapper around the `CloudSetup.launch_instance` function, which is called in parallel using the `multiprocessing` package for efficiency purposes.
4. Creates connections to all the instances using the `connect_instances` function. This is a thin wrapper around the `CloudSetup.connect_instance` function. The function returns an SSH client that can be used to connect to each of the instances.
5. Sets up all the instances using the `setup_instances` function (this step is skipped if the `nocreate` command-line argument is provided). This copies `.boto`, [Cloud/INSTALL.py](INSTALL.py), and [Cloud/INSTALL.sh](INSTALL.sh) from your local `Cloud` folder to the instance. It then runs the `INSTALL.py` script, which runs in an infinite loop, running `INSTALL.sh` with a 6-minute timeout until the `MQLib` executable is built.
6. Splits up the graphs between each machine using the `split_graphs` function. This function first reads in the graphs to be processed. In the case of a `METRICS` run, graphs are read from `Cloud/metric_files.txt` if provided and otherwise from the `mqlibinstances` S3 bucket. In the case of a `FULL` run, graphs are read from [data/runtimes.csv](../data/runtimes.csv). Graphs are then split up evenly between nodes.
7. Starts a run on all EC2 nodes using the `dispatch_and_run` function (this step is skipped if the `nodispatch` command-line argument is provided). Communicates with each EC2 node what graphs it should run by copying a list of graphs to `MQLib/Cloud/GRAPH_FILE` on the node. For `FULL` runs, `GRAPH_FILE` further contains a runtime limit for each instance, and the script further copies [data/heuristics.txt](../data/heuristics.txt) to location `MQLib/Cloud/HEUR_FILE` on each node. Finally, `dispatch_and_run` executes the [Cloud/MQLibMaster.py](MQLibMaster.py) script on each node from the `MQLib/Cloud` folder, which calls the [Cloud/MQLibRunner.py](MQLibRunner.py) script in an infinite loop (see below for details about this script).
8. Terminates once work has been dispatched to all EC2 nodes.

The [Cloud/MQLibRunner.py](MQLibRunner.py) script actually performs all metric computations (for `METRICS` runs) or heuristic executions (for `FULL` runs) on the EC2 nodes. This script is run from the `MQLib/Cloud` folder on each EC2 node and performs a number of steps:

1. Loads all data from the `GRAPH_FILE`, which tells us what graphs to run in what order (and runtime limits in the case of a `FULL` run).
2. In the case of a `FULL` run, loads all data from the `HEUR_FILE`, which lists all heuristics we'll be testing.
3. Connects to the SimpleDB database (`mqlib-metrics` for a `METRICS` run and `mqlib-domain` for a `FULL` run) using `CloudSetup.setup_sdb_domain` and opens log file `Cloud/PROGRESS`.
4. For each graph, uses function `get_graph` to acquire the graph from the `mqlibinstances` S3 bucket and unzip it, storing it in file `Cloud/curgraph`. For a `METRICS` run, computes the metrics if they don't already appear in the `mqlib-metrics` database. For a `FULL` run, executes each heuristic (with fixed seed 144 and the specified runtime limit) if the heuristic has not previously been tested on the graph in the `mqlib-domain` database. Adds the results to the appropriate SimpleDB, using the graph name as the key for the `METRICS` run and the graph name, heuristic code, and repetition number (0) as the key for the `FULL` run. Stored results include the graph name, timestamp, and output of the `MQLib` executable, as well as the heuristic name and repetition id (0) in the case of a `FULL` run.
5. Terminates the EC2 node once all graphs have been processed.

[Dunning et al. (2015)]: http://www.optimization-online.org/DB_FILE/2015/05/4895.pdf
