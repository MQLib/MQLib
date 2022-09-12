# `data` folder contents

The [data/SearchResults_Filtered.txt](SearchResults_Filtered.txt) file contains the search results for the systematic review included in [the paper](../paper/SECM_final.pdf).

The [data/heuristics.txt](heuristics.txt) file contains the set of all heuristics that have been fully tested in the zipped heuristic results in the `data` folder. Currently this is the set of 37 heuristics tested in [the paper](../paper/SECM_final.pdf).

The [data/instance_header_info.csv](instance_header_info.csv) file contains the headers of all the Max-Cut instances in the public `mqlibinstances` S3 bucket. This information can be used to understand the types of instances being tested, as well as the parameters used to generate random instances.

The [data/instances.txt](instances.txt) file contains the file names for the 3,396 problem instances that have been fully tested in the [data/results.zip](results.zip) and [data/results_imgseg.zip](results_imgseg.zip) files and used in [the paper](../paper/SECM_final.pdf).

The [data/metrics.csv](metrics.csv) file contains the metrics returned by the `bin/MQLib` executable for all Max-Cut instances in the public `mqlibinstances` S3 bucket. This file was created by the [Cloud/grabGraphInfo.py](../Cloud/grabGraphInfo.py) script.

The [data/scaling.csv](scaling.csv) file contains the results of the scaling analysis, in which heuristics are run on problem instances of different sizes and densities and their memory usage is measured. This file was created by the [scripts/scaling.py](scripts/scaling.py) script.

The [data/standard.csv](standard.csv) file contains the file names and sources for all problem instances deemed to be in the standard instance library, as described in [the paper](../paper/SECM_final.pdf).

The `data` folder contains the zipped results of computational testing of heuristics, as extracted by the [Cloud/grabFullRuns.py](../Cloud/grabFullRuns.py) script. The following files are currently included:

* The [data/results.zip](results.zip) file contains the best overall solution for the 37 heuristics tested in [the paper](../paper/SECM_final.pdf), across the 3,296 problem instances considered in that paper (that are not for image segmentation). One row is included for each heuristic / instance / seed tuple.
* The [data/results_imgseg.zip](results_imgseg.zip) file contains all best solutions during the run of the 38 heuristics in [the paper](../paper/SECM_final.pdf) across the 100 image segmentation problem instances considered in that work; each heuristic is run for 5 different random seeds. The heuristics were the 37 considered throughout the paper, plus the DESOUSA2013 heuristic.
