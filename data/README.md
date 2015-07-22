# `data` folder contents

The [data/heuristics.txt](heuristics.txt) file contains the set of all heuristics that have been fully tested in the zipped heuristic results in the `data` folder. Currently this is the set of 37 heuristics tested in [Dunning et al. (2015)].

The [data/runtimes.csv](runtimes.csv) file contains the file names and runtime limits for all problem instances that have been fully tested in the zipped heuristics results in the `data` folder, and the [data/metrics.csv](metrics.csv) file contains the metrics returned by the `bin/MQLib` executable for these problem instances. Both of these files are created by the [Cloud/grabGraphInfo.py](../Cloud/grabGraphInfo.py) script. Currently these files describe the 3,398 problem instances tested in [Dunning et al. (2015)] and available in the `mqlibinstances` S3 bucket.

The [data/standard.csv](standard.csv) file contains the file names and sources for all problem instances deemed to be in the standard instance library, as described in [Dunning et al. (2015)].

The `data` folder contains the zipped results of computational testing of heuristics. These files specify the runtimes and solution values of all new best solutions encountered in a set of heuristic runs, as extracted by the [Cloud/grabFullRuns.py](../Cloud/grabFullRuns.py) script. The following files are currently included:

* [data/results.zip](results.zip): All new best solutions for the 37 heuristics tested in [Dunning et al. (2015)](http://www.optimization-online.org/DB_FILE/2015/05/4895.pdf), across the 3,398 problem instances considered in that paper. Additionally, the results of DESOUSA2013 are included for the 100 image segmentation problems.


[Dunning et al. (2015)]: http://www.optimization-online.org/DB_FILE/2015/05/4895.pdf