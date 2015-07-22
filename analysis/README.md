# Analysis folder

This folder contains three R scripts for analyzing the output of evaluations.

* `analyze.R` is the primary script. It generates the tables and figures
  shown in the original paper and can also be used to perform new analysis
  on results generated for new heuristics, instances, and metrics.
* `runtime.R` is the script used to determine the run time limits.
* `storeRF.R` is a helper script used by `analyze.R` to store the trained
  hyperheuristic's random forest models.

Finally `papersummary.csv` is a summary of the attributes of each of papers
that were implemented for the original paper. It is used by `analyze.R` to
generate some tables and figures.
