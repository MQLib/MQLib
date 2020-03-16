################################################################################
# MQLib - analyze.R
# MIT License
# Analyzes results of heuristic evaluation, and produces the results and
# summaries found in the original paper.
################################################################################


################################################################################
# Libraries
#-------------------------------------------------------------------------------
# The analysis uses multiple R packages to aid in the analysis. They can be
# installed by running the following lines:
if (!require("e1071")) {  install.packages("e1071") ; library(e1071)}
if (!require("dplyr")) {  install.packages("dplyr") ; library(dplyr)}  # Aids in manipulating data
if (!require("data.table")) {  install.packages("data.table") ; library(data.table)}  # For faster loading of data
if (!require("ggplot2")) {  install.packages("ggplot2") ; library(ggplot2) }  # Plotting library
if (!require("grid")) {  install.packages("grid") ; library(grid) }  # Used with ggplot2
if (!require("scales")) {  install.packages("scales") ; library(scales) } # Used with ggplot2
if (!require("xtable")) {  install.packages("xtable") ; library(xtable) }  # Produces tables in LaTeX format
if (!require("rpart")) {  install.packages("rpart") ; library(rpart) }  # Produces classification and regression trees
if (!require("rpart.plot")) {  install.packages("rpart.plot") ; library(rpart.plot) }  # Plots those trees
if (!require("randomForest")) {  install.packages("randomForest") ; library(randomForest) }  # Produces random forests
if (!require("mvtnorm")) {  install.packages("mvtnorm") ; library(mvtnorm) }  # Multivariate normal distribution
if (!require("caret")) {  install.packages("caret") ; library(caret) }  # For cross-validating machine learning models
if (!require("FNN")) {  install.packages("FNN") ; library(FNN) }  # efficient k-nearest neighbors
if (!require("tidyr")) {  install.packages("tidyr") ; library(tidyr) }  # For "gather"
if (!require("lazyeval")) {  install.packages("lazyeval") ; library(lazyeval) }  # For "interp"
#-------------------------------------------------------------------------------


################################################################################
# Loading data
#-------------------------------------------------------------------------------
# We will assume we have the following four inputs before proceeding with
# the analysis:
#-------------------------------------------------------------------------------
# Metrics
# Primary metrics data is provided in the data folder (../data/metrics.csv).
# This file was produced for the metrics and instances considered in the
# original paper. The first column is the graph name, and columns 2:59 consist
# of the metrics themselves. Columns 60:71 are the times required to calculate
# the metrics.
initial.metrics <- read.csv("../data/metrics.csv",stringsAsFactors=F)
metrics.times <- data.frame(graphname=initial.metrics$graphname,
                            time=rowSums(initial.metrics[,60:71]))
metrics.all <- initial.metrics[,1:59]  # Just the metrics
# NEW METRICS: if you have new metrics calculated for these original instances
# you should uncomment and modify the following lines.
# new.metrics <- read.csv("../data/newmetrics.csv",stringsAsFactors=F)
# metrics.all <- left_join(metrics.all, new.metrics, by="graphname")
# NEW GRAPHS: if you have new instances and calculated the same metrics on them,
# you should uncomment and modify the following lines to combine them with the
# included metrics.
# new.metrics <- read.csv("../data/newmetrics.csv",stringsAsFactors=F)[,1:59]
# metrics.all <- rbind(metrics.all, new.metrics)
# BOTH NEW GRAPHS AND NEW METRICS: The easiest approach may be simply
# recalculate the metrics on the old graphs and produce a single CSV file.
# 
# Metrics cleanup: if the reported first eigenvalue of the laplacian is smaller than
# the second then flip them since they should be in order.
ev1 <- metrics.all$log_norm_ev1
ev2 <- metrics.all$log_norm_ev2
metrics.all$log_norm_ev1 <- pmax(ev1, ev2)
metrics.all$log_norm_ev2 <- pmin(ev1, ev2)
metrics.all$log_ev_ratio <- abs(metrics.all$log_ev_ratio)
#-------------------------------------------------------------------------------
# Instance headers
# Load the header information for each instance on the server. This information
# can be generated with the `scripts/downloadGraphHeaders.py` script.
initial.instance.headers <- read.csv("../data/instance_header_info.csv",
                                     stringsAsFactors=FALSE)
instance.headers.all <- initial.instance.headers
# NEW GRAPHS: if you have new instances and have created a spreadsheet
# containing the name of the instance (`fname`), the number of nodes (`n`),
# the number of edges (`m`) and the comments from the header of the
# instance (`comments`), then you should uncomment and modify the following
# lines to combine in this new information.
# new.instance.headers <- read.csv("../data/newheaders.csv", stringsAsFactors=FALSE)
# instance.headers.all <- rbind(initial.instance.headers, new.instance.headers)
#-------------------------------------------------------------------------------
# Standard instance list
# Lists all the instances that are part of the "standard instance library" as
# described in the paper. Contains columns `graphname`, `source`, and `problem`
# (currently either `MAX-CUT` or `QUBO`). This file is provided with the source,
# and only needs to be updated if new "standard" instances are added (perhaps
# for a new problem type).
standard <- read.csv("../data/standard.csv",stringsAsFactors=F)
#-------------------------------------------------------------------------------
# Heuristics list
# A list of all heuristics implemented for the original paper and their
# characteristic features.
initial.heuristics <- read.csv("heuristics.csv",stringsAsFactors=F)
heuristics <- initial.heuristics
# NEW HEURISTICS: either combine with the existing file, or uncomment and modify
# the following lines to combine them.
# new.heuristics <- read.csv("newheuristics.csv",stringsAsFactors=F)
# heuristics <- rbind(heuristics, new.heuristics)
#-------------------------------------------------------------------------------
# Evaluation results
# The file ../data/results.zip contains all the results for the 37 heuristics
# tested in the original paper, for 3,296 instances.
# It has columns `timestamp`, `graphname`, `heuristic`, `seed`, `limit`,
# and `objective`, with one row per heuristic run containing the best
# solution found within the runtime limit.

# The file ../data/results_imgseg.zip contains all the results for the
# 100 image segmentation problems. These results have all the columns above
# except `seed` (only one run per instance was performed for each heuristic),
# and it also contains column `runtime`. Each row represents a new best
# solution for the heuristic on the problem instance, so the data taken as
# a whole represents the history of the current best solution for each run.

raw.results <- read.csv(unz("../data/results.zip", "results.csv"), stringsAsFactors = FALSE)
raw.imgseg.results <- read.csv(unz("../data/results_imgseg.zip", "results_imgseg.csv"), stringsAsFactors = FALSE)
# NEW RESULTS: if you have new results for new instances, uncomment and modify
# the following lines to combine them (similarly for image segmentation results).
# new.results <- read.csv("newresults.csv", stringsAsFactors=FALSE)
# raw.results <- rbind(raw.results, new.results)

# Don't group together runs by seed, but compute the deviation for
# each heuristic.
results.byseed <- raw.results %>%
  group_by(graphname) %>%
  mutate(dev = 1 - objective/max(objective))

# Summarize the performance across all replicates for each
# heuristic/instance pairing.
results.all <- raw.results %>%
  group_by(graphname, seed) %>% 
  mutate(rep.rank = rank(-objective, ties.method="min")) %>%
  group_by(graphname, heuristic) %>%
  summarize(best = max(objective),
            average = mean(objective),
            worst = min(objective),
            avgrankrep = mean(rep.rank),
            runs = n())

# If a heuristic was not run for the full number of replicates on,
# an instance, remove the results for that heuristic/instance pair.
replicates <- max(results.all$runs)
if (!all(results.all$runs == replicates)) {
  print(paste("WARNING: There is at least one case where a heuristic",
              "was not run for all", replicates, "replicates on an instance"))
  print(paste("         Removing", sum(results.all$runs != replicates),
              "heuristic/instance pairs"))
  results.all <- results.all %>% filter(runs == replicates)
}

# Compute deviation from the best solution encountered for each instance
results.all <- results.all %>%
  group_by(graphname) %>%
  mutate(rank = rank(-average, ties.method="min"),
         rank2 = rank(-average),
         overallBest = max(best),
         bestDev = ifelse(overallBest == 0, 0, (overallBest - best) / overallBest),
         averageDev = ifelse(overallBest == 0, 0, (overallBest - average) / overallBest),
         worstDev = ifelse(overallBest == 0, 0, (overallBest - worst) / overallBest),
         achievedBest = as.numeric(best == overallBest),
         numHeuristic = n()) %>%
  ungroup() %>%
  select(-overallBest)

# If not all heuristics were run on an instance, remove the instance
num.heuristic <- length(unique(results.all$heuristic))
if (!all(results.all$numHeuristic == num.heuristic)) {
  to.remove <- unique(results.all$graphname[results.all$numHeuristic != num.heuristic])
  print(paste("WARNING: There is at least one case where not all",
              num.heuristic, "heuristics were run on an instance"))
  print(paste("         Removing", length(to.remove), "instances"))
  results.all <- results.all %>% filter(!graphname %in% to.remove)
}

# Some basic summary statistics
print("Section 1: counts of prior heuristics")
print(paste("Number of total instances tested:",
            length(unique(results.all$graphname))))
print(paste("Number of heuristics tested:", num.heuristic))

# For image segmentation, compute the best solution through 1%, 10%, and 100% of
# the runtime limit. Runtime limit < 0 indicates there was no output from the heuristic,
# so these are coded as solution 0.
imgseg.results.all <- raw.imgseg.results %>%
  group_by(graphname, heuristic, seed) %>%
  summarize(obj1 = ifelse(first(limit) < 0, 0, max(objective[runtime <= limit * 0.01])),
            obj10 = ifelse(first(limit) < 0, 0, max(objective[runtime <= limit * 0.1])),
            obj100 = max(objective)) %>%
  group_by(graphname) %>%
  mutate(dev1 = 1 - obj1 / max(obj100),
         dev10 = 1 - obj10 / max(obj100),
         dev100 = 1 - obj100 / max(obj100)) %>%
  ungroup()

################################################################################
# Section 2 and Appenix B: Metric Coverage
#-------------------------------------------------------------------------------
# For input vector `x` (normalized to be in [0, 1]), returns the proportion of
# the line segment [0, 1] that is within `eps` of a vector element.
coverage <- function(x, eps) {
  x <- sort(x)
  lows <- c(x - eps, 1)
  highs <- c(0, x + eps)
  return(1 - sum(pmax(lows - highs, 0)))
}

# Standardized non-binary graph metrics for all instances and just for
# the standard test bed (in both cases limiting to instances with 500+ nodes).
coverage.metrics <- metrics.all[,sapply(metrics.all, function(x) length(unique(x)) > 2)] %>%
  filter(graphname %in% results.all$graphname) %>%
  filter(log_n >= log(500)) %>%
  mutate_if(is.numeric, funs((. - min(.)) / (max(.) - min(.))))
std.coverage.metrics <- coverage.metrics %>%
  filter(graphname %in% standard$graphname)

# Coverage metric for the standard test bed
std.coverage <- std.coverage.metrics %>%
  summarize_if(is.numeric, funs(coverage(., 0.05))) %>%
  gather(metric, coverage)
print("Section 2 and Appendix B coverage information")
print(paste("Standard test bed coverage", mean(std.coverage$coverage),
            "instances", nrow(std.coverage.metrics), "metrics", nrow(std.coverage)))

# Coverage metric for the full expanded test bed
full.coverage <- coverage.metrics %>%
  summarize_if(is.numeric, funs(coverage(., 0.05))) %>%
  gather(metric, coverage)
print(paste("Our full coverage", mean(full.coverage$coverage),
            "instances", nrow(coverage.metrics)))

# Coverage metrics for 1000 iterates of downsampling the expanded library
# to the size of the standard library
set.seed(144)
downsampled <- do.call(rbind, lapply(1:1000, function(idx) {
  coverage.metrics %>%
    sample_n(nrow(std.coverage.metrics)) %>%
    summarize_if(is.numeric, funs(coverage(., 0.05))) %>%
    mutate(iter = idx) %>%
    gather(metric, coverage, -iter)
}))
coverage.per.sample <- downsampled %>%
  group_by(iter) %>%
  summarize(meanCoverage = mean(coverage))
print(paste("Downsampled mean", mean(coverage.per.sample$meanCoverage),
            "2.5% quantile", quantile(coverage.per.sample$meanCoverage, 0.025),
            "97.5% quantile", quantile(coverage.per.sample$meanCoverage, 0.975)))
df.plot <- std.coverage %>%
  inner_join(full.coverage, by="metric", suffix=c("standard", "full")) %>%
  mutate(metric = factor(metric, levels=metric[order(-coveragefull, -coveragestandard)]))
print("Outputting Appendix Figure 1 to coverage.pdf...")
pdf("coverage.pdf", 6, 6)
print(ggplot(df.plot, aes(x=metric, y=coveragefull)) +
  geom_bar(stat="identity") +
    geom_bar(aes(y=coveragestandard), stat="identity", fill="gray") +
    coord_flip() +
    ylab("Coverage") +
    xlab("") +
    theme_bw(base_size=8))
dev.off()

################################################################################
# Appendix A: Summaries of standand instance library
#-------------------------------------------------------------------------------
std.metrics <- left_join(standard, metrics.all, by="graphname")
print("Appendix Section 1 summary of standard instance library:")
print(as.data.frame(std.metrics %>%
  group_by(source) %>%
  summarize(problem = problem[1],
            count = n(),
            min_n = round(min(exp(log_n))),
            max_n = round(max(exp(log_n))),
            min_density = round(min(deg_mean*100), 1),
            max_density = round(max(deg_mean*100), 1))))
#-------------------------------------------------------------------------------


################################################################################
# Figure 1: Overview of how heuristics have been compared historically
# This table and plot do not use data used anywhere else in the figures,
# but are provided here for completeness.
#-------------------------------------------------------------------------------
paper.data <- read.csv("papersummary.csv", stringsAsFactors=F)
print("Section 1 summary information of previous papers")
print(paste("Num Papers:", nrow(paper.data)))
print(paste("Num Papers since 2010:", sum(paper.data$Year >= 2010)))
print(paste("Proportion with publicly accessible source code:", mean(grepl("YES", paper.data$Publish.code))))
print(paste("Don't report processor:", mean(paper.data$Processor == "")))
print(paste("Don't report compiler flags:", mean(paper.data[paper.data$Compiler.Flags != "n/a","Compiler.Flags"] == "")))
paper.data$same.term <- grepl("Same", paper.data$Termination)
paper.data$same.comp <- grepl("Same computer", paper.data$Type)
does.comp <- paper.data[paper.data$CompH > 0,]
print(paste("Compare against another proportion", nrow(does.comp)/nrow(paper.data)))
print(paste("Proportion of comparisons with different termination criteria:", mean(!does.comp$same.term)))
print("Proportion in testing a standard library computed by hand from Graphs column of papersummary.csv")
print(paste("Largest number of publicly accessible instances tested on:", max(paper.data$NumAccessible)))
print(paste("Maximum number compared against", max(paper.data$CompH)))
print("Appendix table 2 constructed by hand from papersummary.csv")
print("Original list of 810 search results can be found in ../data/SearchResults_Filtered.txt")

print("Outputting Figure 1 to broad_eval.pdf...")
p <- ggplot(paper.data, aes(x=CompH)) +
  geom_histogram(binwidth=.5, origin=.25) +
  scale_x_continuous(limits=c(.75, 10.25), breaks=1:10) +
  theme_bw(base_size=9) +
  xlab("Number of Previous Heuristics Compared Against") + ylab("Frequency")
ggsave(filename="broad_eval.pdf", plot=p, width=4, height=3, units="in")
p
#-------------------------------------------------------------------------------


################################################################################
# Figure 2: Scatter plots showing distributions of standard instance libraries
#-------------------------------------------------------------------------------
fig2.data   <- left_join(standard, metrics.all, by="graphname")
fig2.maxcut <- fig2.data %>% filter(problem=="MAX-CUT")
fig2.qubo   <- fig2.data %>% filter(problem=="QUBO")
#-------------------------------------------------------------------------------
# Max-cut plot
print("Outputting Figure 2 (left) to std_maxcut.pdf...")
p <- ggplot(fig2.maxcut, aes(x=exp(log_n), y=deg_mean)) +
  geom_point(size=2) +
  theme_bw(base_size=9) +
  scale_x_log10("Number of nodes",breaks=c(10,100,1000,10000), limits=c(5,30000)) +
  scale_y_continuous("Density",limits=c(0,1))
ggsave(filename="std_maxcut.pdf", plot=p, width=3, height=3, units="in")
p
#-------------------------------------------------------------------------------
print("Outputting Figure 2 (right) to std_qubo.pdf...")
p <- ggplot(fig2.qubo, aes(x=exp(log_n), y=deg_mean)) +
  geom_point(size=2) +
  theme_bw(base_size=9) +
  scale_x_log10("Number of variables",breaks=c(10,100,1000,10000), limits=c(5,30000)) +
  scale_y_continuous("Density",limits=c(0,1))
ggsave(filename="std_qubo.pdf", plot=p, width=3, height=3, units="in")
p
#-------------------------------------------------------------------------------

################################################################################
# Figure 4: Sources of instances in the expanded instance library
#-------------------------------------------------------------------------------
headers.keep <- instance.headers.all %>%
  filter(fname %in% results.all$graphname)
if (nrow(headers.keep) != length(unique(results.all$graphname))) {
  print("Header information is incomplete")
} else {
  # If you have expanded the set of instances to other source types, then 
  # you should expand this dictionary, which indicates what the search for in
  # the instance comments for each instance source.
  identifiers <- data.frame(search.string = c("G-set problem instance", "Beasley instances from ORLIB", "Culberson graph. ", "random graph using networkx library in python", "COMMENT: biqmac_all/mac/rudy data set", "COMMENT: biqmac_all/biq/be/ data set", "STP File, STP", "random graph using rudy library in python. ", "\\.col", "Section Comment, Name", "COMMENT: biqmac_all/mac/ising data set", "SECTION Comment, Name", "instances generated by helmberg and rendl", "SECTION Comments, Name", "this set contains 30 instances described in festa et al", "Glover, Kochenberger, and Alidaee \\(1998\\) ", "Pardalos and Rogers \\(1989\\) ", "srgraphs\\.html", "SECTION Comment Name", "SECTION Comment , Name", "funkybee\\.narod\\.ru", "www\\.dharwadker\\.org", "Processed from file bqp", "COMMENT: biqmac_all/biq/beasley data set", "From Dimacs 7", "TYPE: TSP", "COMMENT: c FILE: ", "www\\.proin\\.ktu\\.lt", "\\(GKA\\) instances from ORLIB", "33D32945 GRP File", "hiv-2"),
                            inst.source = c("G-set", "ORLIB", "Culberson", "networkx", "mac/rudy", "biq/be", "steinlib", "rudy", "GraphColoring", "steinlib", "mac/ising", "steinlib", "H+R 2000", "steinlib", "F 2002", "GKA", "PR 1998", "Internet", "steinlib", "steinlib", "Internet", "Internet", "ORLIB", "biq/beasley", "Dimacs 7", "TSPLIB", "GraphColoring", "Palubeckis", "GKA", "steinlib", "steinlib"),
                            source.type = c("Library", "Library", "Random", "Random", "Library", "Library", "Library", "Random", "Library", "Library", "Library", "Library", "Library", "Library", "Library", "Library", "Library", "Library", "Library", "Library", "Library", "Library", "Library", "Library", "Library", "Library", "Library", "Library", "Library", "Library", "Library"))
  headers.keep$sourceNum <- NA
  for (row in seq_len(nrow(identifiers))) {
    headers.keep$sourceNum[grep(identifiers$search.string[row],
                                headers.keep$comments)] <- row
  }
  if (any(is.na(headers.keep$sourceNum))) {
    print("Not all instances could be identified based on comments")
  } else {
    headers.keep$inst.source <- identifiers$inst.source[headers.keep$sourceNum]
    headers.keep$source.type <- identifiers$source.type[headers.keep$sourceNum]
    print("Source of our instances (see Appendix Section 1) -- see data/instance_header_info.csv for detailed parameters")
    print(table(headers.keep$inst.source))
    print("High-level type of our instances (leftmost two numbers in Figure 4):")
    print(table(headers.keep$source.type))
  }
}

################################################################################
# Figure 5: Scatter plots showing distribution of expanded instance library
#-------------------------------------------------------------------------------
fig5.data <- metrics.all %>%
  left_join(standard, by="graphname") %>%
  mutate(type=ifelse(is.na(problem), "Added", problem)) %>%
  arrange(type)  # So added instances are plotted below Maxcut and QUBO
#-------------------------------------------------------------------------------
print("Outputting Figure 5 (left) to stdvsexp_EV.pdf...")
p <- ggplot(fig5.data, aes(x=log_norm_ev2, y=log_ev_ratio, color=type)) +
  geom_point(size=2) +
  scale_color_manual(guide=FALSE, values=c("darkgrey","blue","red")) +
  theme_bw(base_size=9) +
  scale_x_continuous("log(normalized 2nd Laplacian eigenvalue)") +
  scale_y_continuous("log(ratio of top 2 Laplacian eigenvalues)")
ggsave(filename="stdvsexp_EV.pdf", plot=p, width=3, height=3, units="in")
p
#-------------------------------------------------------------------------------
print("Outputting Figure 5 (left) to stdvsexp_nvsdense.pdf...")
p <- ggplot(fig5.data, aes(x=exp(log_n), y=deg_mean, color=type)) +
  geom_point(size=2) +
  scale_color_manual(guide=FALSE, values=c("darkgrey","blue","red")) +
  theme_bw(base_size=9) +
  scale_x_log10("Number of nodes",breaks=c(10,100,1000,10000)) +
  scale_y_continuous("Density")
ggsave(filename="stdvsexp_nvsdense.pdf", plot=p, width=3, height=3, units="in")
p
#-------------------------------------------------------------------------------


################################################################################
# Identification of "interesting" instances
# "Interesting" is defined to be the instances for which less than half, i.e.
# floor(37/2)=18 (if using the original 37) or fewer heuristics matched the
# best solution across all their replicates. This is equivalent to saying the
# median of the 37 heuristics' average deviations is greater than zero.
int.graphs <- (results.all %>%
  group_by(graphname) %>%
  summarize(propPerfect = mean(worstDev == 0)) %>%
  filter(propPerfect <= 0.5))$graphname

# The number of interesting instances is 2735 for original paper's instance
# library and 37 heuristics
print("Number of interesting instances (Section 4.2):")
print(length(int.graphs))

# Collect the results for only interesting instances
int.results <- results.all %>% filter(graphname %in% int.graphs)

# A quick CART tree to see if we can understand what sorts of graphs were removed
understand.int.graphs <- results.all %>%
  group_by(graphname) %>%
  summarize(uninteresting = as.numeric(mean(worstDev == 0) > 0.5)) %>%
  inner_join(metrics.all, by="graphname")
print("Breakdown of small vs. large in terms of interestingness (first paragraph of Section 4.2):")
print(understand.int.graphs %>%
  group_by(exp(log_n) >= 100) %>%
  summarize(propUninteresting = mean(uninteresting)))
larger.understand <- understand.int.graphs %>% filter(exp(log_n) >= 100)
uninteresting.mod <- rpart(uninteresting~.-graphname, data=larger.understand,
                           minbucket=200)
prp(uninteresting.mod, varlen = 100)
#-------------------------------------------------------------------------------
# Section 4: times and costs
# Compute the number of runs for each instance at each runtime limit
num.runs.rtlim <- raw.results %>%
  group_by(graphname, limit) %>%
  summarize(runs = n())
# Compute CPU-hours of runs (remembering that for each instance/runtime
# limit pair we also did a baseline run)
cpu.hours <- with(num.runs.rtlim, sum(limit*(runs+1))) / 3600
print("Section 4.1: cost of evaluating heuristics")
print(paste("Total runtime (CPU-years)", cpu.hours / 24 / 365.25,
            "Per heuristic (CPU-days)", cpu.hours / num.heuristic / 24,
            "Average total time/node (days)", cpu.hours / 24 / 60,
            "Average single heuristic time/node (hours)", cpu.hours / num.heuristic / 60,
            "Total cost ($)", cpu.hours * 0.067,
            "Cost/heuristic ($)", cpu.hours / num.heuristic * 0.067))

#-------------------------------------------------------------------------------
# Section 7: times and costs
# Total human hours to run (per heuristic), assuming 20 nodes
print("Section 7 costs to test a new heuristic:")
print("Wall-clock hours, limiting to interesting instances")
print(sum(raw.results$limit[raw.results$graphname %in% int.graphs]) / num.heuristic / 3600 / 20)
print("Wall-clock hours, across all instances")
print(sum(raw.results$limit) / num.heuristic / 3600 / 20)
# Total human dollars to run (per heuristic, $0.07 per hour)
print("Cost, limiting to interesting instances")
print(sum(raw.results$limit[raw.results$graphname %in% int.graphs]) / num.heuristic / 3600 * 0.07)
print("Cost, across all instances")
print(sum(raw.results$limit) / num.heuristic / 3600 * 0.07)
#-------------------------------------------------------------------------------


################################################################################
# Compare best-known to optimal in cases where we know optimal
#-------------------------------------------------------------------------------
known.opt <- results.all %>%
  group_by(graphname) %>%
  summarize(bestHeur = max(best)) %>%
  inner_join(standard, by="graphname") %>%
  filter(!is.na(optimal))
print("Section 4.2: comparison to optimal when it's known")
print(paste(sum(known.opt$bestHeur == known.opt$optimal),
            "instances had an optimal best heuristic solution and",
            sum(known.opt$bestHeur != known.opt$optimal), "didn't"))
print("Section 4.2: Summary of instance props for standard test bed with known opt:")
print("  Remember converting QUBO->MAX-CUT adds 1 to the size")
print(table(round(exp(metrics.all$log_n))[metrics.all$graphname %in% known.opt$graphname]))
#-------------------------------------------------------------------------------


################################################################################
# Generate summaries of results (Figure 7, Table 2)
#-------------------------------------------------------------------------------
exp.summary <- int.results %>%
  group_by(heuristic) %>%
  summarize(firstequal  = mean(rank == 1) * 100,
            firststrict = mean(rank2 == 1) * 100,
            pctbest     = mean(achievedBest) * 100,
            worstdev    = mean(worstDev) * 100,
            meandev     = mean(averageDev) * 100,
            bestdev     = mean(bestDev) * 100,
            avgrankrep  = mean(avgrankrep)) %>%
  arrange(desc(firstequal), desc(firststrict), meandev)
exp.pretty.summary <- exp.summary %>%
  mutate(heuristic = heuristics$abbreviation[match(heuristic, heuristics$heuristic)])
#-------------------------------------------------------------------------------
# Table 1: Summary of results
print("Table 1 latex output:")
summary.latex <- xtable(exp.pretty.summary)
digits(summary.latex) <- 1
print(summary.latex, floating=TRUE, include.rownames=FALSE)
#-------------------------------------------------------------------------------
# Figure 7: FIRST-EQUAL (breaks and limits have been hand-selected and may need
#           to be adjusted based on new results)
merge.final <- int.results %>%
  filter(graphname %in% standard$graphname) %>%
  group_by(heuristic) %>%
  summarize(firstequal  = mean(rank == 1) * 100,
            firststrict = mean(rank2 == 1) * 100,
            meandev     = mean(averageDev) * 100,
            avgrankrep  = mean(avgrankrep)) %>%
  left_join(exp.summary, by="heuristic", suffix=c(".std", ".exp")) %>%
  select(heuristic, firstequal.std, firstequal.exp, firststrict.std,
         firststrict.exp, meandev.std, meandev.exp, avgrankrep.std,
         avgrankrep.exp) %>%
  mutate_if(is.numeric, funs("rank"=rank(-.))) %>%
  arrange(desc(firstequal.std))
print("Outputting Figure 7 (top-left) to firstequal_stdvsexp.pdf...")
p <-
  ggplot(merge.final, aes(x=firstequal.std,y=firstequal.exp)) +
  geom_abline(slope=1,color="gray") + scale_linetype_identity() +
  geom_point(size=2) +
  scale_x_continuous(name="% first-equal, standard instance lib.", breaks=seq( 0,40,10), limits=c(0,40)) +
  scale_y_continuous(name="% first-equal, expanded instance lib.", breaks=seq( 0,40,10), limits=c(0,40)) +
  theme_bw(base_size=9)
ggsave(filename="firstequal_stdvsexp.pdf", plot=p, width=3, height=3, units="in")
p
#-------------------------------------------------------------------------------
# Figure 7: FIRST-STRICT (breaks and limits have been hand-selected and may need
#           to be adjusted based on new results)
print("Outputting Figure 7 (top-right) to firststrict_stdvsexp.pdf...")
p = ggplot(merge.final, aes(x=firststrict.std, y=firststrict.exp)) + 
  geom_abline(slope=1,color="gray") + scale_linetype_identity() +
  geom_point(size=2) +
  scale_x_continuous(name="% first-strict, standard instance lib.", breaks=seq( 0,30,5),limits=c(0,32)) +
  scale_y_continuous(name="% first-strict, expanded instance lib.", breaks=seq( 0,30,5),limits=c(0,32)) +
  theme_bw(base_size=9)
ggsave(filename="firststrict_stdvsexp.pdf", plot=p, width=3, height=3, units="in")
p
#-------------------------------------------------------------------------------
# Figure 7: MEAN-DEV (breaks and limits have been hand-selected and may need
#           to be adjusted based on new results)
print("Outputting Figure 7 (bottom-left) to meandev_stdvsexp.pdf...")
p = ggplot(merge.final, aes(x=meandev.std, y=meandev.exp)) + 
  geom_abline(slope=1,color="gray") + scale_linetype_identity() +
  geom_point(size=2) +
  scale_x_log10(name="% mean deviation, standard instance lib.",
                limits=c(0.1,100), breaks=c(.1, 1, 10, 100)) +
  # , labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(name="% mean deviation, expanded instance lib.",
                limits=c(0.1,100), breaks=c(.1, 1, 10, 100)) +
  # , labels = trans_format("log10", math_format(10^.x))) +
  theme_bw(base_size=9)
ggsave(filename="meandev_stdvsexp.pdf", plot=p, width=3, height=3, units="in")
p
#-------------------------------------------------------------------------------
# Figure 7: MEAN-RANK (breaks and limits have been hand-selected and may need
#           to be adjusted based on new results)
print("Outputting Figure 7 (bottom-right) to avgrank_stdvsexp.pdf...")
p = ggplot(merge.final, aes(x=avgrankrep.std, y=avgrankrep.exp)) + 
  geom_abline(slope=1,color="gray") + scale_linetype_identity() +
  geom_point(size=2) +
  scale_x_continuous(name="Average rank, standard instance lib.", limits=c(1, 37)) +
  # , labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(name="Average rank, expanded instance lib.", limits=c(1, 37)) +
  # , labels = trans_format("log10", math_format(10^.x))) +
  theme_bw(base_size=9)
ggsave(filename="avgrank_stdvsexp.pdf", plot=p, width=3, height=3, units="in")
p
#-------------------------------------------------------------------------------
print("Section 4.2: MER99LS summary on standard and new instances")
mer99ls <- subset(merge.final, heuristic=="MERZ1999GLS")
print(paste("new instances first strict percentage",
            mer99ls$firststrict.exp, "rank", mer99ls$firststrict.exp_rank))
print(paste("standard instances first strict percentage",
            mer99ls$firststrict.std, "first equal percentage",
            mer99ls$firstequal.std, "first equal rank",
            sum(merge.final$firstequal.std > mer99ls$firstequal.std)+1))
#-------------------------------------------------------------------------------
print("Section 4.2: first-strict changes between standard and expanded libraries")
print(paste("0 first-strict on standard but positive on expanded:",
            sum(merge.final$firststrict.std == 0 & merge.final$firststrict.exp > 0),
            "; >= 5%", sum(merge.final$firststrict.std == 0 & merge.final$firststrict.exp >= 5)))

################################################################################
# Section 6: Generate and evaluate hyperheuristic
# As described in the paper, we use one random forest model for each heuristic.
#-------------------------------------------------------------------------------

# We won't build a random forest model for a heuristic if it ties the best
# solution on 10 or fewer instances. This prevents issues with rare outcomes
# when fitting the random forests and corresponds to thinking of the
# algorithm portfolio as always assigning probability 0 that these
# heuristics will match the best, which is not a terrible approximation.
hh.removed.heuristics <- (int.results %>%
  group_by(heuristic) %>%
  summarize(numBest = sum(rank == 1)) %>%
  filter(numBest < 10))$heuristic
print(paste("Removing", paste(hh.removed.heuristics, collapse=", "),
            "from the hyperheuristic"))
results.for.hh <- int.results %>%
  filter(!heuristic %in% hh.removed.heuristics)
# Set a random seed to help with reproducibility
set.seed(123)
# Identify the set of all interesting graphs, then split them into two sets
graphnames <- sort(unique(int.graphs))
train.graphs <- sample(graphnames, 0.7*length(graphnames))
test.graphs <- graphnames[!graphnames %in% train.graphs]
print(paste("Section 6.1: Train size", length(train.graphs), "test size",
            length(graphnames)-length(train.graphs)))
# Remove metrics we don't want to use. We currently use all of them.
#to.remove <- c("...")
#red.metrics <- metrics.all[,!(names(metrics) %in% to.remove)]
red.metrics <- metrics.all
#-------------------------------------------------------------------------------
# Build a random forest model for each heuristic, and evaluate it on the test
# set. Will also generate and save the random forest models, which can then
# be moved to /hhdata/ if desired.
# Load the random forest saving functionality
if (!dir.exists("fittedTrees")) {
  dir.create("fittedTrees")
}
source("storeRF.R")

# Set this to true if you want to load cached versions instead of re-running the
# cross-validation
load.cached <- TRUE

# Iterate over every heuristic
hh.raw.result.list <- lapply(split(results.for.hh, results.for.hh$heuristic), function(df){
  print(df$heuristic[1])
  
  # Which rows of the per-heuristic dataframe correspond to training graphs?
  train.rows <- which(df$graphname %in% train.graphs)
  
  # INDEPENDENT VARIABLES
  # Get the metrics for all graphs in the dataframe
  X <- red.metrics[match(df$graphname, red.metrics$graphname),]
  # Remove the non-numeric column graphname
  X$graphname <- NULL
  # Convert to a matrix for machine learning
  X <- as.matrix(X)
  # Training testing split
  X.train  <- X[ train.rows,]
  X.test   <- X[-train.rows,]
  X.imgseg <- red.metrics %>% filter(graphname %in% imgseg.results.all$graphname)
  imgseg.graphs <- X.imgseg$graphname
  X.imgseg <- X.imgseg %>% select(-graphname)
  
  # DEPENDENT VARIABLE
  # 1 iff heuristic is first equal on this graph, 0 otherwise
  y <- as.factor(as.numeric(df$rank == 1))
  # Training testing split
  y.train <- y[ train.rows]
  y.test  <- y[-train.rows]
  
  # RANDOM FOREST
  set.seed(144)
  if (load.cached && file.exists(paste0("fittedTrees/", df$heuristic[1],".rds"))) {
    cvResults <- NULL
    finalModel <- readRDS(paste0("fittedTrees/", df$heuristic[1],".rds"))
  } else {
    rf.train <- train(X.train, y.train, method="rf",
                      trControl=trainControl(method="cv",number=5))
    finalModel <- rf.train$finalModel
    cvResults <- rf.train$results
    # Save the best forest
    store.rf(finalModel, paste0("fittedTrees/", df$heuristic[1],".rf"))
    saveRDS(finalModel, paste0("fittedTrees/", df$heuristic[1],".rds"))
  }
  
  # Store training-set predictions
  pred.train <- data.frame(heuristic = df$heuristic[1],
                           graphname = df$graphname[train.rows],
                           pred      = predict(finalModel, newdata=X.train, "prob"),
                           y         = y.train)
  
  # EVALUATION ON TEST SET
  pred <- predict(finalModel, newdata=X.test, "prob")
  if (!is.null(cvResults)) {
    print(paste("Acc train  rf:", max(rf.train$results$Accuracy)))
  }
  print(paste("Acc test base:", max(table(y.test))/length(y.test)))
  print(paste("Acc test   rf:", confusionMatrix(predict(finalModel, newdata=X.test), y.test)$overall[1]))
  pred.test <- data.frame(heuristic = df$heuristic[1],
                          graphname = df$graphname[-train.rows],
                          pred      = pred,  # Actually pred.0 and pred.1
                          y         = y.test)
  
  # Evaluation on image segmentation instances
  pred.imgseg <- data.frame(graphname = imgseg.graphs,
                            pred      = predict(finalModel, newdata=X.imgseg,
                                                type="prob")[,"1"],
                            heuristic = df$heuristic[1])
  
  # Variable importance
  imp <- importance(finalModel)
  imp <- data.frame(variable         = row.names(imp),
                    MeanDecreaseGini = imp,
                    heuristic        = df$heuristic[1])

  # Return a list of results
  return(list(pred.train  = pred.train,
              pred.test   = pred.test,
              pred.imgseg = pred.imgseg,
              imp         = imp))
})
print("Random forest models stored in fittedTrees/*.rf and *.rds")

# Combine test set results for each heuristic
hh.raw.result <- do.call(rbind, lapply(hh.raw.result.list, "[[", "pred.test"))
# Cache results
print("Random forest predictions for each instance in results_hh.csv")
write.csv(hh.raw.result, file="results_hh.csv", row.names=F)

# Combine training set results for each heuristic
hh.raw.result.train <- do.call(rbind, lapply(hh.raw.result.list, "[[", "pred.train"))

# Combine image segmentation results for each heuristic
hh.raw.imgseg <- do.call(rbind, lapply(hh.raw.result.list, "[[", "pred.imgseg"))

# Combine the variable importance scores from the random forest models
hh.var.imp <- do.call(rbind, lapply(hh.raw.result.list, "[[", "imp"))

#-------------------------------------------------------------------------------
# Summarize variable importance measures across the fitted trees
print("Outputting appendix table 3 (variable importance summary) latex to var_imp_summary.txt...")
var.imp.summary <- hh.var.imp %>%
  group_by(heuristic) %>%
  mutate(MDGrank = rank(-MeanDecreaseGini)) %>%
  ungroup() %>%
  group_by(variable) %>%
  summarize(meanMDG = round(mean(MeanDecreaseGini), 1),
            meanMDGrank = round(mean(MDGrank), 1),
            pctTop10 = round(mean(MDGrank <= 10) * 100)) %>%
  arrange(meanMDGrank) %>%
  select(variable, meanMDGrank, meanMDG, pctTop10)
var.imp.latex <- xtable(var.imp.summary)
print(var.imp.latex, floating=TRUE, include.rownames=FALSE,
      file="var_imp_summary.txt")
print(paste("Section 6.1: Number of features not in the top 10 for any random forest:", sum(var.imp.summary$pctTop10 == 0)))





#-------------------------------------------------------------------------------
# Given a set of instances and all heuristics to run on those instances, returns
# the expected mean deviation and expected probability of a best solution for each.
analyze.k.heur <- function(data, results = results.byseed) {
  # Add a unique ID for each row
  data <- data %>% mutate(idx = seq_len(nrow(data)))
  
  # Compute the probability of each deviation if we only ran using the first
  # heuristic
  condensed <- data %>%
    inner_join(results, by=c("graphname", "heur1"="heuristic")) %>%
    group_by(graphname, idx, dev) %>%
    summarize(prob = n() / replicates) %>%
    ungroup()
  # One-by-one merge in information about the next heuristics to predict the
  # probability of each deviation being the minimum one encountered.
  hnum <- 2
  while(paste0("heur", hnum) %in% names(data)) {
    col <- paste0("heur", hnum)
    this.by <- setNames(c("graphname", "heuristic"), c("graphname", col))
    condensed <- condensed %>%
      inner_join(data %>% select(idx, match(col, names(.))), by="idx") %>%
      inner_join(results, by=this.by) %>%
      mutate(dev = pmin(dev.x, dev.y)) %>%
      group_by(graphname, idx, dev) %>%
      summarize(prob = sum(prob) / replicates)
    hnum <- hnum + 1
  }
  summarized <- condensed %>%
    group_by(graphname, idx) %>%
    summarize(expectedDev = sum(dev * prob),
              probBest = sum((dev == 0) * prob))
  data %>%
    inner_join(summarized, by=c("graphname", "idx")) %>%
    select(-idx)
}

# The naive approach orders by first-equal; analyze those results
naive.ordering <- (exp.summary %>%
  arrange(desc(firstequal)))$heuristic
naive.info <- do.call(rbind, lapply(1:8, function(k) {
  naive.choice <- data.frame(graphname = test.graphs,
                             stringsAsFactors=FALSE)
  for (i in 1:k) {
    naive.choice[,paste0("heur", i)] <- naive.ordering[i]
  }
  naive.results <- analyze.k.heur(naive.choice)
  naive.results$k <- k
  naive.results %>%
    select(graphname, k, expectedDev, probBest)
}))
naive.summary <- naive.info %>%
  group_by(k) %>%
  summarize(expectedDev = mean(expectedDev) * 100,
            pctBest = mean(probBest) * 100) %>%
  mutate(type = "Naive")

print("Section 6.1 summary of performance of a method that orders results")
print("  by first-equal rate and takes the top N")
print(naive.summary)

# The naive hyperheuristic approach orders by the hyperheuristic's
# prediction of being first-equal.
set.seed(144)
hh.ranked <- hh.raw.result %>%
  mutate(heuristic = as.character(heuristic),
         graphname = as.character(graphname)) %>%
  group_by(graphname) %>%
  mutate(hhprobrank=rank(pred.0,ties.method="random")) %>%
  ungroup()
print(paste("Number heuristics used >= 1 time in 1-heuristic portfolio:",
            length(unique(hh.ranked$heuristic[hh.ranked$hhprobrank == 1]))))
hh.info <- do.call(rbind, lapply(1:8, function(k) {
  hh.choice <- hh.ranked %>%
    filter(hhprobrank <= k) %>%
    mutate(colname = paste0("heur", hhprobrank)) %>%
    select(graphname, heuristic, colname) %>%
    spread(colname, heuristic)
  hh.results <- analyze.k.heur(hh.choice)
  hh.results$k <- k
  hh.results %>%
    select(graphname, k, expectedDev, probBest)
}))
hh.summary <- hh.info %>%
  group_by(k) %>%
  summarize(expectedDev = mean(expectedDev) * 100,
            pctBest = mean(probBest) * 100) %>%
  mutate(type = "Random forest (independent)")

print("Hyperheuristic summary (select N best heuristics according to RF models):")
print(hh.summary)

# For each k-tuple, compute the minimum, independent, and maximum
# probability of first-equal for each passed instance using the
# hyperheuristic passed probabilities as well as whether one of
# the heuristics was actually first-equal for that instance.
hh.summarize <- function(hh.dat, k) {
  hh.base <- hh.dat %>%
    mutate(graphname = as.character(graphname),
           heur = as.character(heuristic),
           outcome = as.numeric(as.character(y)),
           low = pred.1, mid = pred.1, high = pred.1) %>%
    select(graphname, heur, outcome, low, mid, high)
  hh.combined <- hh.base %>%
    rename(heur1 = heur)
  for (i in 2:k) {
    col <- paste0("heur", i)
    prevCol <- paste0("heur", i-1)
    hh.combined <- hh.combined %>%
      inner_join(hh.base, by="graphname")
    # A few base R operations for readability
    hh.combined[,col] <- hh.combined$heur
    hh.combined <- hh.combined[hh.combined[,prevCol] < hh.combined[,col],]
    # Finish off the chain by combining in the new information
    hh.combined <- hh.combined %>%
      mutate(outcome = pmax(outcome.x, outcome.y),
             low = pmax(low.x, low.y),
             mid = 1 - (1 - mid.x) * (1 - mid.y),
             high = high.x + high.y) %>%
      select(-outcome.x, -outcome.y, -low.x, -low.y, -mid.x, -mid.y,
             -high.x, -high.y, -heur)
  }
  hh.combined
}

kbest.list <- lapply(2:3, function(k) {
  # For each tuple, compute the average outcome and average low, mid, and
  # high assessment across the training set
  train.summary <- hh.summarize(hh.raw.result.train, k) %>%
    group_by_(.dots = paste0("heur", 1:k)) %>%
    summarize(meanOutcome = mean(outcome),
              meanLow = mean(low),
              meanMid = mean(mid),
              meanHigh = mean(high))
  
  # For each tuple, compute "alpha", which takes value -1 if meanOutcome is
  # at or below meanLow, takes value 0 at meanMid (linear interpolation in between),
  # and takes value 1 at or abovemeanHigh (linear interpolation in between).
  train.summary <- train.summary %>%
    mutate(low.val = pmax(pmin((meanOutcome - meanLow) / (meanMid - meanLow), 1), 0),
           high.val = pmax(pmin((meanOutcome - meanMid) / (meanHigh - meanMid), 1), 0)) %>%
    mutate(alpha = -1 + low.val + high.val) %>%
    select(-low.val, -high.val)
  
  # For the test set, compute the low/mid/high for each topic and instance and then
  # merge in alpha and use it to compute the estimated probability of success for the tuple
  test.combined <- hh.summarize(hh.raw.result, k) %>%
    inner_join(train.summary %>% select(-meanOutcome, -meanLow, -meanMid, -meanHigh),
               by=paste0("heur", 1:k)) %>%
    mutate(pred = low + pmin(alpha+1, 1) * (mid - low) + pmax(alpha, 0) * (high - mid))
  
  # Select the heuristic to run for each instance as the one with the highest prediction,
  # breaking ties randomly
  set.seed(144)
  kbest.choice <- test.combined %>%
    group_by(graphname) %>%
    mutate(rank = rank(-pred, ties.method="random")) %>%
    filter(rank == 1) %>%
    ungroup()

  # Determine how well the selected tuples did on the test set
  kbest.results <- analyze.k.heur(kbest.choice)
  kbest.results$k <- k
  kbest.results <- kbest.results %>%
    select(graphname, k, expectedDev, probBest)

  # Return both `train.summary`, which provides `alpha` for each tuple, as well as
  # `kbest.results`, which shows how we actually did on the test set.
  list(tuple.info = train.summary,
       results    = kbest.results)
})
# The k=1 case is equivalent to the normal hyperheuristic, so merge in those
# results from hh.info
kbest.results <- rbind(hh.info %>% filter(k == 1), 
                       do.call(rbind, lapply(kbest.list, "[[", "results")))
kbest.summary <- kbest.results %>%
  group_by(k) %>%
  summarize(expectedDev = mean(expectedDev) * 100,
            pctBest = mean(probBest) * 100) %>%
  mutate(type = "Random forest (correlated)")

#-------------------------------------------------------------------------------
# Generate Figure 10
top.N.results <- rbind(naive.summary, hh.summary)
print("Outputting Figure 10 (left) to hh_vs_base_top8_fe.pdf...")
p <- ggplot(top.N.results, aes(x=k, y=pctBest, linetype=type)) +
  geom_line() +
  geom_point() +
  scale_linetype_discrete(guide=FALSE) +
  scale_x_continuous("Number of heuristics selected",breaks=seq(1,8),minor_breaks=1:8) +
  scale_y_continuous("Mean Chance of Returning Best Solution (%)") +
  theme_bw(base_size=9)
ggsave(filename="hh_vs_base_top8_fe.pdf", plot=p, width=3, height=3, units="in")
p
print("Outputting Figure 10 (right) to hh_vs_base_top8_dev.pdf...")
p <- ggplot(top.N.results, aes(x=k, y=expectedDev, linetype=type)) +
  geom_line() +
  geom_point() +
  scale_linetype_discrete(guide=FALSE) +
  scale_x_continuous("Number of heuristics selected",breaks=seq(1,8),minor_breaks=1:8) +
  scale_y_continuous("Mean Expected Deviation (%)") +
  theme_bw(base_size=9)
ggsave(filename="hh_vs_base_top8_dev.pdf", plot=p, width=3, height=3, units="in")
p
# Now do a version that has labels and that also includes the k-best approach
top.N.results.expanded <- rbind(naive.summary, hh.summary, kbest.summary)
print("Outputting an attempt at a smarter hyperheuristic (not in paper) to hh_vs_base_top8_fe_expanded.pdf and hh_vs_base_top8_dev_expanded.pdf...")
p <- ggplot(top.N.results.expanded, aes(x=k, y=pctBest, linetype=type, col=type)) +
  geom_line() +
  geom_point() +
  scale_x_continuous("Number of heuristics selected",breaks=seq(1,8),minor_breaks=1:8) +
  scale_y_continuous("Percentage First Equal") +
  theme_bw(base_size=9)
ggsave(filename="hh_vs_base_top8_fe_expanded.pdf", plot=p, width=6, height=3, units="in")
p
p <- ggplot(top.N.results.expanded, aes(x=k, y=expectedDev, linetype=type, col=type)) +
  geom_line() +
  geom_point() +
  scale_x_continuous("Number of heuristics selected",breaks=seq(1,8),minor_breaks=1:8) +
  scale_y_continuous("Mean deviation (%)") +
  theme_bw(base_size=9)
ggsave(filename="hh_vs_base_top8_dev_expanded.pdf", plot=p, width=6, height=3, units="in")
p

#-------------------------------------------------------------------------------
# Perform the out-of-sample analysis of the image segmentation application.
set.seed(144)
imgseg.ranked <- hh.raw.imgseg %>%
  mutate(heuristic = as.character(heuristic),
         graphname = as.character(graphname)) %>%
  group_by(graphname) %>%
  mutate(hhprobrank = rank(-pred, ties.method="random")) %>%
  ungroup()
imgseg.hh.choice <- imgseg.ranked %>%
  filter(hhprobrank == 1) %>%
  mutate(colname = "heur1") %>%
  select(graphname, heuristic, colname) %>%
  spread(colname, heuristic)
imgseg.results.100 <- imgseg.results.all %>%
  mutate(dev = dev100)
imgseg.hh.summary.100 <- analyze.k.heur(imgseg.hh.choice, imgseg.results.100)
imgseg.results.10 <- imgseg.results.all %>%
  mutate(dev = dev10)
imgseg.hh.summary.10 <- analyze.k.heur(imgseg.hh.choice, imgseg.results.10)
print("Section 6.2 image segmentation results:")
print(paste("HH dev at 100% RTlim", mean(imgseg.hh.summary.100$expectedDev) * 100,
            "at 10% RTlim", mean(imgseg.hh.summary.10$expectedDev) * 100))
desousa.choice <- data.frame(graphname = imgseg.hh.choice$graphname,
                             heur1 = "DESOUSA2013",
                             stringsAsFactors = FALSE)
desousa.summary.100 <- analyze.k.heur(desousa.choice, imgseg.results.100)
desousa.summary.10 <- analyze.k.heur(desousa.choice, imgseg.results.10)
print(paste("desousa dev at 100% RTlim", mean(desousa.summary.100$expectedDev) * 100,
            "at 10% RTlim", mean(desousa.summary.10$expectedDev) * 100))
print("Section 6.2: Summary of metric computation runtimes for image segmentation:")
print(summary(metrics.times$time[grepl("imgseg", metrics.times$graphname)]))

hh.info <- do.call(rbind, lapply(1:8, function(k) {
  hh.choice <- hh.ranked %>%
    filter(hhprobrank <= k) %>%
    mutate(colname = paste0("heur", hhprobrank)) %>%
    select(graphname, heuristic, colname) %>%
    spread(colname, heuristic)
  hh.results <- analyze.k.heur(hh.choice)
  hh.results$k <- k
  hh.results %>%
    select(graphname, k, expectedDev, probBest)
}))
hh.summary <- hh.info %>%
  group_by(k) %>%
  summarize(expectedDev = mean(expectedDev) * 100,
            pctBest = mean(probBest) * 100) %>%
  mutate(type = "Random forest (independent)")

print("Hyperheuristic summary:")
print(hh.summary)

EXCLUDE.IMGSEG <- FALSE
if (EXCLUDE.IMGSEG) {
  hh.ranked.imgseg <- hh.ranked %>%
    filter(grepl("imgseg",hh.ranked$graphname)) %>%
    filter(hhprobrank<=1)
  table(hh.ranked.imgseg$heuristic)
  summary(hh.ranked.imgseg$rank)
  summary(hh.ranked.imgseg$dev)
  
  hh.ranked.imgseg <- hh.ranked %>%
    filter(grepl("imgseg",hh.ranked$graphname)) %>%
    filter(heuristic=="PALUBECKIS2004bMST2")
  summary(hh.ranked.imgseg$rank)
  summary(hh.ranked.imgseg$dev)
}
#-------------------------------------------------------------------------------


################################################################################
# Section 5.1: Understanding a single heuristic
# Here we build trees and plots for all heuristics as shown in the paper.
#-------------------------------------------------------------------------------
# For each heuristic, compute a simple CART tree predicting instance
# difficulty, returning summary information.
heuristic.names <- sort(unique(int.results$heuristic))
if (!dir.exists("treeOutput")) {
  dir.create("treeOutput")
}
print("Outputting trees and plots of top 2 splits to treeOutput folder...")
print("  Figure 8 (top) is treeOutput/tree_LAGUNA2009CE.pdf and treeOutput/plot_LAGUNA2009CE.pdf")
print("  Figure 8 (middle) is treeOutput/tree_DUARTE2005.pdf and treeOutput/plot_DUARTE2005.pdf")
print("  Figure 8 (bottom) is treeOutput/tree_MERZ2002GREEDYKOPT.pdf and treeOutput/plot_MERZ2002GREEDYKOPT.pdf")
tree.info <- do.call(rbind, lapply(heuristic.names, function(heur) {
  # Fit a CART tree predicting the normalized deviation using all metrics
  res.sng <- int.results %>%
    filter(heuristic == heur) %>%
    select(graphname, avgrankrep) %>%
    left_join(red.metrics, by="graphname")
  mod <- rpart(avgrankrep ~ . - graphname, data=res.sng, minbucket=100,
               maxdepth = 3)

  # In this section we will be extracting information from CART trees 
  # returned by rpart (`mod` above) to determine the top 2 most important
  # variables and to partition the instance space defined by these two
  # variables.
  # 
  # `mod$frame` contains information about each node in the tree; if
  # it is a split node then `var` will be set to the split variable, and
  # otherwise it is set to "<leaf>". `ncompete+nsurrogate` is the number
  # of competing and surrogate splits at that split point. The row number
  # indicates the position in the tree: the root is 1, its left and right
  # children are 2 and 3, the next row is 4, 5, 6, 7, etc. The parent of
  # node i is node floor(i/2), the left child of node i is 2*i and the
  # right child of node i is node 2*i+1.
  # 
  # `mod$splits` contains information about all the splits, in addition to
  # the competing and surrogate splits. It is in the same order as
  # `mod$frame`, so the split at the root node appears first, followed by
  # all its competing and surrogate splits, then the next split appears
  # followed by all its competing and surrogate splits, and so on.
  # 
  # First, we identify the first two unique variables to have been used as a
  # split in the tree.
  breaks <- as.character(mod$frame$var)
  breaks <- unique(breaks[breaks != "<leaf>"])
  if (length(breaks) >= 2) {
    var1 <- min(breaks[1:2])
    var2 <- max(breaks[1:2])
  } else {
    var1 <- ifelse(length(breaks) >= 1, breaks[1], NA)
    var2 <- NA
  }
  
  # Tree R^2
  SSE <- sum((predict(mod) - res.sng$avgrankrep)^2)
  SST <- sum((mean(res.sng$avgrankrep) - res.sng$avgrankrep)^2)

  # Tree (manually drawn in paper)
  pdf(paste0("treeOutput/tree_", heur, ".pdf"), width=4, height=3)
  prp(mod, digits=4, varlen=0, type=3, fallen.leaves=TRUE, tweak=0.8,
      node.fun=function(x, labs, digits, varlen) {
        labs <- lapply(labs, function(l) { substr(l,1,4) })
        return(paste(labs, "\nn=",x$frame$n, sep=""))
      })
  dev.off()
  
  if (is.na(var1) || is.na(var2)) {
    return(data.frame(heuristic = heur, r2 = 1 - SSE/SST, var1, var2,
                      stringsAsFactors=FALSE))
  }
  
  # Make the node number and corresponding row in mod$splits easily
  # accessible in mod$frame
  mod$frame$node <- as.numeric(row.names(mod$frame))
  mod$frame$split.row <-
    head(cumsum(c(1, mod$frame$ncompete + mod$frame$nsurrogate +
                     (mod$frame$var != "<leaf>"))), -1)
    
  # Compute the splits that we will plot (anything left after we
  # remove any split on a variable other than the two we keep, also
  # removing the subtree below any such splits)
  to.rem <- mod$frame$node[!mod$frame$var %in% c(var1, var2)]
  for (iter in 1:3) {
    to.rem <- unique(c(to.rem, 2*to.rem, 2*to.rem+1))
  }
  frame.rows.keep <- !mod$frame$node %in% to.rem
  split.info <- mod$frame[frame.rows.keep,]
  
  # Add information about the split cutoff value and direction for the
  # left side of each split.
  split.info$cutoff <- mod$splits[split.info$split.row, "index"]
  split.info$left.dir <- ifelse(mod$splits[split.info$split.row, "ncat"] < 0, "<", ">=")
  
  # Compute the segments to be plotted for the CART model
  mins <- c(min(res.sng[,var1]), min(res.sng[,var2]))
  names(mins) <- c(var1, var2)
  maxes <- c(max(res.sng[,var1]), max(res.sng[,var2]))
  names(maxes) <- c(var1, var2)

  # Each row in split.info specifies a split on a variable; this function returns
  # the ranges of the other variable for which the split in row "row" of split.info
  # is valid by looking at the ancestor splits in the CART tree that split the other
  # variable as well as the range of the raw data for the other variable. "mins" and
  # "maxes" are named vectors of the min and max of both variables of interest in
  # the full dataset.
  bounds <- function(row, split.info, mins, maxes) {
    # To find all ancestors, keep dividing by 2 and rounding down
    # until you reach 1
    reduction <- Reduce("%/%", rep(2, 10), split.info$node[row], accumulate=TRUE)
    reduction <- reduction[reduction >= 1]
    
    # You are on the left of a parent if your node number is two times theirs
    ancestor.nodes = tail(reduction, -1)
    is.left = head(reduction, -1) == 2*tail(reduction, -1)
    
    # Grab some information about the ancestors from split.info
    ancestor.var <- split.info$var[match(ancestor.nodes, split.info$node)]
    ancestor.cutoff <- split.info$cutoff[match(ancestor.nodes, split.info$node)]
    ancestor.left.dir <- split.info$left.dir[match(ancestor.nodes, split.info$node)]
    ancestor.dir <- ifelse(is.left, ancestor.left.dir,
                           ifelse(ancestor.left.dir == "<", ">=", "<"))
    
    # You are lower-bounded by ancestor splits on the other variable that are ">="
    # for your side and upper-bounded by ancestor splits on the other variable that
    # are "<" for your side
    this.var <- split.info$var[row]
    lower <- max(c(mins[names(mins) != this.var],
                   ancestor.cutoff[ancestor.var != this.var & ancestor.dir == ">="]))
    upper <- min(c(maxes[names(maxes) != this.var],
                   ancestor.cutoff[ancestor.var != this.var & ancestor.dir == "<"]))
    c("lower" = lower, "upper" = upper)
  }
  
  # Apply the "bounds" function to each split of interest to get the
  # lower and upper bounds at which the split applies.
  split.info$lower <- sapply(seq_len(nrow(split.info)), function(row) {
    bounds(row, split.info, mins, maxes)["lower"]
  })
  split.info$upper <- sapply(seq_len(nrow(split.info)), function(row) {
    bounds(row, split.info, mins, maxes)["upper"]
  })

  p <- ggplot(res.sng, aes_string(x=var1,y=var2,color="avgrankrep")) +
         scale_color_continuous(guide=FALSE, low="blue", high="red",
                                limits=range(int.results$avgrankrep)) +
         geom_point(size=1) +
         theme_bw(base_size=9)
  for (row in seq_len(nrow(split.info))) {
    if (split.info$var[row] == var1) {
      p <- p + geom_segment(x = split.info$cutoff[row],
                            y = split.info$lower[row],
                            xend = split.info$cutoff[row],
                            yend = split.info$upper[row],
                            col = "gray", size=0.5)
    } else {
      p <- p + geom_segment(x = split.info$lower[row],
                            y = split.info$cutoff[row],
                            xend = split.info$upper[row],
                            yend = split.info$cutoff[row],
                            col = "gray", size=0.5)
    }
  }

  p
  ggsave(filename=paste0("treeOutput/plot_", heur, ".pdf"),
         plot=p, width=3, height=3, units="in")

  # Return summary information about the tree, including the
  # first two variables used for splits.
  data.frame(heuristic = heur, r2 = 1 - SSE/SST, var1, var2,
             stringsAsFactors=FALSE)
}))

# Summary information of the R^2 of trees
print("Section 5.1: Summary information of trees' R^2 values:")
print(summary(tree.info$r2))

# Summary information about heuristics that share the same two most
# important variables.
tree.info %>%
  group_by(var1, var2) %>%
  filter(n() >= 2)

# Latex table of heuristic names and the R^2 of the CART tree
tree.output <- tree.info %>%
  mutate(heuristic = heuristics$abbreviation[match(heuristic, heuristics$heuristic)],
         r2 = round(r2, 2)) %>%
  arrange(heuristic) %>%
  select(heuristic, r2) %>%
  mutate(figure = paste0("Figure", heuristic))
tree.summary.latex <- xtable(tree.output)
# digits(summary.latex) <- 1
print("Outputting Appendix Table 4 to tree_summary.txt...")
print(tree.summary.latex, floating=TRUE, include.rownames=FALSE,
      sanitize.text.function = NULL, file="tree_summary.txt")
#-------------------------------------------------------------------------------


################################################################################
# Section 5.2: Understanding heuristic ideas
# Here we provide a generic function to create the density plots as shown
# in the paper. By changing the heuristic features data, or by adding new
# heuristics/instances, the results can be updated
#-------------------------------------------------------------------------------
# See what labels look like they might make sense
res.class <- results.all %>%
  filter(rank2 == 1.0,
         graphname %in% int.graphs) %>%
  left_join(heuristics, by="heuristic")
print(paste("Section 5.2 number of instances with exactly one best heuristic", nrow(res.class)))
# The desireable property is that one feature is not dominant over the others
print("Breakdown by various heuristic properties:")
print("Classification:")
print(table(res.class$classification))  # Iterated Local Search
print("Hybrid:")
print(table(res.class$Hybrid))  # Balanced
print("Initialization:")
print(table(res.class$Initialization))  # Balanced
print("Population:")
print(table(res.class$Population))  # Balanced
print("Evolutionary:")
print(table(res.class$Evolutionary))  # No vs 1+2Parent
print("Memory:")
print(table(res.class$Memory))  # Balanced
print("Perturbation:")
print(table(res.class$Perturbation))  # Balanced
print("Problem:")
print(table(res.class$problem))  # Not very balanced
#-------------------------------------------------------------------------------
# Define the plotting function
nnplot <- function(xdata, ydata, outcome, pts, filename, k, rad, xlab, ylab, title) {
  df <- data.frame(x=xdata, y=ydata, outcomes=outcome)
  
  # Determine area
  xmin <- ifelse(min(df$x) < 0, min(df$x) * 1.01, min(df$x) * 0.99)
  xmax <- ifelse(max(df$x) < 0, max(df$x) * 0.99, max(df$x) * 1.01)
  ymin <- ifelse(min(df$y) < 0, min(df$y) * 1.01, min(df$y) * 0.99)
  ymax <- ifelse(max(df$y) < 0, max(df$y) * 0.99, max(df$y) * 1.01)

  # Setup grid
  xvals <- seq(xmin, xmax, length.out=pts)
  yvals <- seq(ymin, ymax, length.out=pts)
  points <- as.matrix(expand.grid(x=xvals, y=yvals))

  # Compute the mean outcome among the k nearest neighbors of each plot point
  # as well as the number of these k neighbors that are within "rad" euclidean
  # distance of the plot point.
  knn.info <- get.knnx(cbind(df$x, df$y), points, k)
  nn <- rowMeans(matrix(df$outcomes[knn.info$nn.index], nrow(points)))
  n.range <- rowSums(knn.info$nn.dist < rad)
 
  # Desired mapping:
  # 1.0  ->  (1,0,0)
  # 0.5  ->  (0,0,0)
  # 0.0  ->  (0,0,1)
  # Corrected for the number of points in a ball
  nn.red <- ifelse(nn>=0.5,2*(nn-0.5),0) * (n.range/k)
  nn.blue <- ifelse(nn<0.5,2*(0.5-nn),0) * (n.range/k)
  rgba.cols <- rgb(1-nn.red, 1-pmax(nn.blue,nn.red), 1-nn.blue)

  png(filename)
  plot(c(xmin,xmax), range(df$y), type = "n", xlab = xlab, ylab = ylab, main=title, xaxs="i", yaxs="i", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, xaxt="n", yaxt="n")
  # These are the labels used in the paper, and can be uncommented if desired
  # They will not be correct for other metrics
   axis(1, at=c(1,2,3,4,5), labels=c(expression(10^1),
                                     expression(10^2),
                                     expression(10^3),
                                     expression(10^4),
                                     expression(10^5)))
   axis(2, at=c(2,3,4,5,6,7), labels=c(expression(10^2),
                                   expression(10^3),
                                   expression(10^4),
                                   expression(10^5),
                                   expression(10^6),
                                   expression(10^7)))
  rasterImage(t(matrix(rgba.cols, nrow=length(xvals))[,(length(xvals):1)]), xmin, ymin, xmax, ymax)
  points(df$x, df$y, pch=16, cex=0.2, col="black")
  dev.off()
}
#-------------------------------------------------------------------------------
# EVOLUTIONARY
res.EVO <- res.class %>%
  select(graphname,Evolutionary) %>%
  mutate(Evolutionary=ifelse(Evolutionary=="No","No","Yes")) %>%
  mutate(Evolutionary=as.factor(Evolutionary)) %>%
  left_join(red.metrics, by="graphname")
  # %>% filter(weight_mean>=-3.5e-01)
print("Outputting Figure 9 (left) to nn_evo.png...")
nnplot(log10(exp(res.EVO$log_n)), log10(exp(res.EVO$log_m)), as.numeric(res.EVO$Evolutionary=="Yes"),
       400, "nn_evo.png", k=20, rad=1, xlab="Number of Nodes", ylab="Number of Edges", title="Evolutionary Algorithm")
#-------------------------------------------------------------------------------
# MEMORY (TABU SEARCH)
res.MEM <- res.class %>%
  select(graphname,Memory) %>%
  mutate(Memory=as.factor(Memory)) %>%
  left_join(red.metrics, by="graphname")
  # %>% filter(weight_mean>=-3.5e-01)
print("Outputting Figure 9 (right) to nn_mem.png...")
nnplot(log10(exp(res.MEM$log_n)), log10(exp(res.MEM$log_m)), as.numeric(res.MEM$Memory=="Yes"),
       400, "nn_mem.png", k=20, rad=1, xlab="Number of Nodes", ylab="Number of Edges", title="Tabu Search")
#-------------------------------------------------------------------------------

################################################################################
# Appendix Section 3: Scaling analysis
# Here, we plot the scaling analysis for sparse and dense graphs.
#-------------------------------------------------------------------------------
# Load up the memory scaling data
scaling <- read.csv("../data/scaling.csv")

# Plot the scaling for the complete graph case
scaling.complete <- scaling %>%
  filter(graph == "complete.txt") %>%
  mutate(heurPlot = ifelse(heuristic == "PALUBECKIS2004bMST5", "PAL04T5", "Other"))
print("Outputting Appendix Figure 2 (left) to scalingDense.pdf...")
p <- ggplot(scaling.complete, aes(x=size, y=memusg/1024, group=heuristic, color=heurPlot)) +
  geom_point() +
  geom_line() +
  scale_x_log10(name="Number of Nodes", breaks=c(10, 30, 100, 300, 1000, 3000)) +
  scale_y_log10(name="Memory Usage (MB)") +
  scale_color_manual(name="Heuristic", values=c("PAL04T5"="red", "Other"="black")) +
  theme_bw(base_size=9) + 
  theme(legend.position=c(0, 1.043), legend.justification=c(0, 1))
ggsave(filename="scalingDense.pdf", plot=p, width=3, height=3, units="in")
p

# Plot the scaling for the sparse graph case
scaling.sparse <- scaling %>%
  filter(graph == "ER.txt") %>%
  mutate(heurPlot = ifelse(heuristic == "LAGUNA2009CE", "LAG09CE",
                    ifelse(heuristic == "LAGUNA2009HCE", "LAG09HCE",
                    ifelse(heuristic == "PARDALOS2008", "PAR08", "Other"))))
print("Outputting Appendix Figure 2 (right) to scalingSparse.pdf...")
p <- ggplot(scaling.sparse, aes(x=size, y=memusg/1024, group=heuristic, color=heurPlot)) +
  geom_point() +
  geom_line() +
  scale_x_log10(name="Number of Nodes", breaks=c(10, 30, 100, 300, 1000, 3000, 10000, 30000)) +
  scale_y_log10(name="Memory Usage (MB)") +
  scale_color_manual(name="Heuristic", values=c("LAG09CE"="red", "LAG09HCE"="purple",
                                                "PAR08"="blue", "Other"="black")) +
  theme_bw(base_size=9) + 
  theme(legend.position=c(0, 1.043), legend.justification=c(0, 1))
ggsave(filename="scalingSparse.pdf", plot=p, width=3, height=3, units="in")
p




