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
# install.packages(c("dplyr","data.table","ggplot2","grid","scales","xtable"))
# install.packages(c("rpart","rpart.plot","randomForest","mvtnorm","caret"))
library(dplyr)  # Aids in manipulating data
library(data.table)  # For faster loading of data
library(ggplot2)  # Plotting library
library(grid)  # Used with ggplot2
library(scales)  # Used with ggplot2
library(xtable)  # Produces tables in LaTeX format
library(rpart)  # Produces classification and regression trees
library(rpart.plot)  # Plots those trees
library(randomForest)  # Produces random forests
library(mvtnorm)  # Multivariate normal distribution
library(caret)  # For cross-validating machine learning models
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
# the metrics, which we will not use here.
initial.metrics <- read.csv("../data/metrics.csv",stringsAsFactors=F)[,1:59]
metrics.all <- initial.metrics
# NEW METRICS: if you have new metrics calculated for these original instances
# you should uncomment and modify the following lines.
# new.metrics <- read.csv("../data/nrewmetrics.csv",stringsAsFactors=F)
# metrics.all <- left_join(metrics.all, new.metrics, by="graphname")
# NEW GRAPHS: if you have new instances and calculated the same metrics on them,
# you should uncomment and modify the following lines to combine them with the
# included metrics.
# new.metrics <- read.csv("../data/nrewmetrics.csv",stringsAsFactors=F)[,1:59]
# metrics.all <- rbind(metrics.all, new.metrics)
# BOTH NEW GRAPHS AND NEW METRICS: The easiest approach may be simply
# recalculate the metrics on the old graphs and produce a single CSV file.
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
# tested in the original paper, for 3398 instances and 100 image segmentation
# problems (which are treated slightly different).
# It has columns `timestamp`, `graphname`, `heuristic`, `limit`, `objective`,
# and `runtime`, with one row per improving solution found for each heuristic
# and instance. This script doesn't use the raw data: it only needs the final
# improving solution. As a result, the script will cache the reduced results
# instead of computing them every time.
# If you have new results, they will be loaded at the same time. Make sure to
# clear the cached file (results_reduced.csv) otherwise new results will not
# be detected and loaded.
if (file.exists("results_reduced.csv")) {
  print("LOADING CACHED RESULTS (results_reduced.csv)")
  results.all <- read.csv("results_reduced.csv", stringsAsFactors=F)
} else {
  # Load the old results. They are stored in a zip file, so we will dynamically
  # unzip them as we load them. This operation will take a while to run.
  raw.results <- read.csv(unz("../data/results.zip","results.csv"), stringsAsFactors=FALSE)
  # NEW RESULTS: if you have new results for new instances, uncomment and modify
  # the following lines to combine them.
  # new.results <- read.csv("newresults.csv", stringsAsFactors=FALSE)
  # raw.results <- rbind(raw.results, new.results)
  
  # Process results to extract best solution and calculate relative performance
  results.all <- raw.results %>%
    filter(runtime <= limit) %>%
    group_by(graphname,heuristic) %>%
    summarize(limit = limit[1],
              obj = max(objective)) %>%
    group_by(graphname) %>%
    mutate(rank  = rank(-obj, ties.method="min"),
           rank2 = rank(-obj),
           dev   = (max(obj) - obj)/max(obj)) %>%
    ungroup()
  
  # Cache the results
  write.csv(results.all, file="results_reduced.csv", row.names=FALSE)
}
#-------------------------------------------------------------------------------


################################################################################
# Optional: removal of image segmentation graphs
# In the paper, the analysis excluded the image segmentation graphs entirely.
# To recreate the results of the paper, they must be excluded here before
# proceeding. We will use the EXCLUDE.IMGSEG flag later for the hyperheuristic.
#-------------------------------------------------------------------------------
EXCLUDE.IMGSEG <- TRUE
#-------------------------------------------------------------------------------
if (EXCLUDE.IMGSEG) {
  metrics.imgseg <- subset(metrics.all,  grepl("imgseg",metrics.all$graphname))
  metrics        <- subset(metrics.all, !grepl("imgseg",metrics.all$graphname))
  results.imgseg <- subset(results.all,  grepl("imgseg",results.all$graphname))
  results        <- subset(results.all, !grepl("imgseg",results.all$graphname))
} else {
  metrics <- metrics.all
  results <- results.all
}
#-------------------------------------------------------------------------------


################################################################################
# Appendix A: Summaries of standand instance library
#-------------------------------------------------------------------------------
std.metrics <- left_join(standard, metrics, by="graphname")
std.metrics %>%
  group_by(source) %>%
  summarize(problem = problem[1],
            count = n(),
            min_n = min(exp(log_n)),
            max_n = max(exp(log_n)),
            min_density = min(deg_mean*100),
            max_density = max(deg_mean*100))
#-------------------------------------------------------------------------------


################################################################################
# Figure 1: Overview of how heuristics have been compared historically
# This table and plot do not use data used anywhere else in the figures,
# but are provided here for completeness.
#-------------------------------------------------------------------------------
paper.data <- read.csv("papersummary.csv", stringsAsFactors=F)
print(paste("Num Papers:", nrow(paper.data)))
print(paste("Num Papers since 2010:", sum(paper.data$Year >= 2010)))
paper.data$same.term <- grepl("Same", paper.data$Termination)
paper.data$same.comp <- grepl("Same computer", paper.data$Type)
does.comp <- paper.data[paper.data$CompH > 0,]
print("Breakdown of fair comparison (rows are same hardware, cols are same termination criteria)")
table(does.comp$same.comp, does.comp$same.term) / nrow(does.comp)
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
fig2.data   <- left_join(standard, metrics, by="graphname")
fig2.maxcut <- fig2.data %>% filter(problem=="MAX-CUT")
fig2.qubo   <- fig2.data %>% filter(problem=="QUBO")
#-------------------------------------------------------------------------------
# Max-cut plot
p <- ggplot(fig2.maxcut, aes(x=exp(log_n), y=deg_mean)) +
  geom_point(size=2) +
  theme_bw(base_size=9) +
  scale_x_log10("Number of nodes",breaks=c(10,100,1000,10000), limits=c(5,30000)) +
  scale_y_continuous("Density",limits=c(0,1))
ggsave(filename="std_maxcut.pdf", plot=p, width=3, height=3, units="in")
p
#-------------------------------------------------------------------------------
p <- ggplot(fig2.qubo, aes(x=exp(log_n), y=deg_mean)) +
  geom_point(size=2) +
  theme_bw(base_size=9) +
  scale_x_log10("Number of variables",breaks=c(10,100,1000,10000), limits=c(5,30000)) +
  scale_y_continuous("Density",limits=c(0,1))
ggsave(filename="std_qubo.pdf", plot=p, width=3, height=3, units="in")
p
#-------------------------------------------------------------------------------


################################################################################
# Figure 5: Scatter plots showing distribution of expanded instance library
#-------------------------------------------------------------------------------
fig5.data <- metrics %>%
  left_join(standard, by="graphname") %>%
  mutate(type=ifelse(is.na(problem), "Added", problem)) %>%
  arrange(type)  # So added instances are plotted below Maxcut and QUBO
#-------------------------------------------------------------------------------
p <- ggplot(fig5.data, aes(x=weight_stdev, y=assortativity, color=type)) +
  geom_point(size=2) +
  scale_color_manual(guide=FALSE, values=c("darkgrey","blue","red")) +
  theme_bw(base_size=9) +
  scale_x_continuous("Standard Deviation of Weights (normalized)") +
  scale_y_continuous("Degree Assortativtiy")
ggsave(filename="stdvsexp_assvsweight.pdf", plot=p, width=3, height=3, units="in")
p
#-------------------------------------------------------------------------------
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
# best solution. This is equivalent to saying the median deviation is greater
# than zero.
int.graphs <- results %>%
  group_by(graphname) %>%
  mutate(ismax=(obj==max(obj))) %>%
  summarize(interest=sum(ismax)) %>%
  filter(interest<=floor(nrow(heuristics)/2))

# The number of interesting instances is 1968 for original papers instance
# library and 37 heuristics
nrow(int.graphs)

# Collect the results for only interesting instances
int.results <- results %>% filter(graphname %in% int.graphs$graphname)

# Check that median deviation > 0 (following code should return 0)
int.results %>%
  group_by(graphname) %>%
  summarize(meddev=median(dev)) %>%
  mutate(zerodev=meddev==0) %>%
  summarize(sum(zerodev))
#-------------------------------------------------------------------------------
# Section 7: times and costs
# Total human hours to run (per heuristic)
sum(int.results$limit) / nrow(heuristics) / 3600 / 20
sum(results$limit) / nrow(heuristics) / 3600 / 20
# Total human dollars to run (per heuristic, $0.07 per hour)
sum(int.results$limit) / nrow(heuristics) / 3600 * 0.07
sum(results$limit) / nrow(heuristics) / 3600 * 0.07
#-------------------------------------------------------------------------------


################################################################################
# Generate summaries of results (Figure 7, Table 3)
#-------------------------------------------------------------------------------
exp.summary <- int.results %>%
  group_by(heuristic) %>%
  summarize(firstequal  = sum(rank ==1) / n() * 100,
            firststrict = sum(rank2==1) / n() * 100,
            meandev     = mean(dev) * 100) %>%
  arrange(desc(firstequal))

std.summary <- int.results %>%
  filter(graphname %in% standard$graphname) %>%
  group_by(heuristic) %>%
  summarize(firstequal  = sum(rank ==1) / n() * 100,
            firststrict = sum(rank2==1) / n() * 100,
            meandev     = mean(dev) * 100) %>%
  arrange(desc(firstequal))

merge.summary <- inner_join(exp.summary, std.summary, by="heuristic")
merge.final <- data.frame(
  Heuristic   = merge.summary$heuristic,
  FullFirstEq = merge.summary$firstequal.x,
  FullFirstSt = merge.summary$firststrict.x,
  FullMeanDev = merge.summary$meandev.x,
  StdFirstEq  = merge.summary$firstequal.y,
  StdFirstSt  = merge.summary$firststrict.y,
  StdMeanDev  = merge.summary$meandev.y)
merge.final$Heuristic <- as.character(merge.final$Heuristic)
#-------------------------------------------------------------------------------
# Table 3: Summary of results
summary.latex <- xtable(merge.final)
digits(summary.latex) <- 1
print(summary.latex, floating=TRUE, include.rownames=FALSE)
#-------------------------------------------------------------------------------
# Figure 7: FIRST-EQUAL
p <-
  ggplot(merge.final, aes(x=StdFirstEq,y=FullFirstEq)) +
  geom_abline(slope=1,color="gray") + scale_linetype_identity() +
  geom_point(size=2) +
  scale_x_continuous(name="% first-equal, standard instance lib.", breaks=seq( 0,50,10), limits=c(0,50)) +
  scale_y_continuous(name="% first-equal, expanded instance lib.", breaks=seq( 0,50,10), limits=c(0,50)) +
  theme_bw(base_size=9)
ggsave(filename="firstequal_stdvsexp.pdf", plot=p, width=3, height=3, units="in")
p
#-------------------------------------------------------------------------------
# Figure 7: FIRST-STRICT
p = ggplot(merge.final, aes(x=StdFirstSt, y=FullFirstSt)) + 
  geom_abline(slope=1,color="gray") + scale_linetype_identity() +
  geom_point(size=2) +
  scale_x_continuous(name="% first-strict, standard instance lib.", breaks=seq( 0,16,2),limits=c(0,16)) +
  scale_y_continuous(name="% first-strict, expanded instance lib.", breaks=seq( 0,16,2),limits=c(0,16)) +
  theme_bw(base_size=9)
ggsave(filename="firststrict_stdvsexp.pdf", plot=p, width=3, height=3, units="in")
p
#-------------------------------------------------------------------------------
# Figure 7: MEAN-DEV
p = ggplot(merge.final, aes(x=StdMeanDev, y=FullMeanDev)) + 
  geom_abline(slope=1,color="gray") + scale_linetype_identity() +
  geom_point(size=2) +
  scale_x_log10(name="% mean deviation, standard instance lib.",
                limits=c(10^-1.0,20),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(name="% mean deviation, expanded instance lib.",
                limits=c(10^-1.0,20),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw(base_size=9)
ggsave(filename="meandev_stdvsexp.pdf", plot=p, width=3, height=3, units="in")
p
#-------------------------------------------------------------------------------


################################################################################
# Section 6: Generate and evaluate hyperheuristic
# As described in the paper, we use one random forest model for each heuristic.
#-------------------------------------------------------------------------------
# Get the subset of results that we need:
#  - if EXCLUDE.IMGSEG, add it back in
#  - if not, no need for anything special
if (EXCLUDE.IMGSEG) {
  results.for.hh <- results.all %>%
    filter(graphname %in% int.graphs$graphname | grepl("imgseg",graphname),
           heuristic != "DESOUSA2013")
} else {
  results.for.hh <- int.results
}
# Because LAGUNA2009CE does so poorly, models struggle to produce anything
# sensible. We thus exclude it to prevent any issues - we can think of it as
# always predicting that LAGUNA2009CE isn't a good choice.
results.for.hh <- filter(results.for.hh, heuristic != "LAGUNA2009CE")
# Set a random seed to help with reproducibility
set.seed(123)
# Identify the set of all interesting graphs, then split them into two sets
graphnames <- sort(unique(int.graphs$graphname))
train.graphs <- sample(graphnames, 0.7*length(graphnames))
# Remove metrics we don't want to use. We currently use all of them.
#to.remove <- c("...")
#red.metrics <- metrics.all[,!(names(metrics) %in% to.remove)]
red.metrics <- metrics.all
#-------------------------------------------------------------------------------
# Build a random forest model for each heuristic, and evaluate it on the test
# set. Will also generate and save the random forest models, which can then
# be moved to /hhdata/ if desired.
# Load the random forest saving functionality
source("storeRF.R")
# Only runs if no cached results.
if (file.exists("results_hh.csv")) {
  print("LOADING CACHED HYPERHEURISTIC RESULTS (results_hh.csv)")
  hh.raw.result <- read.csv("results_hh.csv",stringsAsFactors=F)
} else {
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
    X.train <- X[ train.rows,]
    X.test  <- X[-train.rows,]
    
    # DEPENDENT VARIABLE
    # 1 iff heuristic is first equal on this graph, 0 otherwise
    y <- as.factor(as.numeric(df$rank == 1))
    # Training testing split
    y.train <- y[ train.rows]
    y.test  <- y[-train.rows]
    
    # RANDOM FOREST
    set.seed(144)
    rf.train <- train(X.train, y.train, method="rf",
                      trControl=trainControl(method="cv",number=2))
    # Save the best forest
    store.rf(rf.train$finalModel, paste(df$heuristic[1],".rf",sep=""))
    
    # EVALUATION
    pred <- predict(rf.train$finalModel, newdata=X.test, "prob")
    print(paste("Acc train  rf:", max(rf.train$results$Accuracy)))
    print(paste("Acc test base:", max(table(y.test))/length(y.test)))
    print(paste("Acc test   rf:", confusionMatrix(predict(rf.train$finalModel, newdata=X.test), y.test)$overall[1]))
    # Return a dataframe of results
    return(data.frame(heuristic = df$heuristic[1],
                      graphname = df$graphname[-train.rows],
                      pred      = pred,  # Actually pred.0 and pred.1
                      y         = y.test))
  })
  # Combine results for each heuristic
  hh.raw.result <- bind_rows(hh.raw.result.list)
  # Cache results
  write.csv(hh.raw.result, file="results_hh.csv", row.names=F)
}
#-------------------------------------------------------------------------------
# Turn probabilities into ranks for each graph, and combine with "true" ranks
hh.ranked <- hh.raw.result %>%
  group_by(graphname) %>%
  arrange(desc(pred.1)) %>%
  mutate(hhprobrank=rank(pred.0,ties.method="random")) %>%
  ungroup() %>%
  left_join(results.all, by=c("heuristic","graphname"))
# Merge in information about the heuristics
# At this point, we remove the image segmentation results if we are trying to
# reproduce the paper results
if (EXCLUDE.IMGSEG) {
  hh.base.ranked <- hh.ranked %>%
    filter(!grepl("imgseg",hh.ranked$graphname)) %>%
    left_join(exp.summary %>% mutate(baserank=rank(-firstequal)), by="heuristic")
} else {
  hh.base.ranked <- hh.ranked %>%
    left_join(exp.summary %>% mutate(baserank=rank(-firstequal)), by="heuristic")
}
# Now for each n=1,.,,8 take the top n for HH, and overall top n (baseline)
top.N.results <- rbind_all(lapply(1:8, function (n) {
  # Baseline
  base.top.n <- hh.base.ranked %>%
    filter(baserank <= n)
  base.top.n.res <- base.top.n %>%
    group_by(graphname) %>%
    summarize(baserank = min(rank),
              basedev  = min(dev),
              baseobj  = max(obj))
  
  # HH
  hh.top.n <- hh.base.ranked %>%
    filter(hhprobrank <= n)
  hh.top.n.res <- hh.top.n %>%
    group_by(graphname) %>%
    summarize(hhrank = min(rank),
              hhdev  = min(dev),
              hhobj  = max(obj))
  
  # Key metrics:
  return(data.frame(
    n=n,
    hhfirstequals   = table(hh.top.n.res$hhrank)[1]/nrow(hh.top.n.res)*100,
    hhdevmean       = mean(hh.top.n.res$hhdev)*100,
    hhdev95         = sort(hh.top.n.res$hhdev)[0.95*nrow(hh.top.n.res)]*100,
    hhnumheurused   = length(table(hh.top.n$heuristic)),
    basefirstequals = table(base.top.n.res$baserank)[1]/nrow(base.top.n.res)*100,
    basedevmean     = mean(base.top.n.res$basedev)*100,
    basedev95       = sort(base.top.n.res$basedev)[0.95*nrow(base.top.n.res)]*100
    ))
}))
#-------------------------------------------------------------------------------
# Generate Figure 10
p = ggplot(top.N.results, aes(x=n)) + 
  geom_line(aes(y=hhfirstequals,linetype="Hyperheuristic")) + geom_point(aes(y=hhfirstequals)) +
  geom_line(aes(y=basefirstequals,linetype="Baseline")) + geom_point(aes(y=basefirstequals)) +
  scale_linetype_discrete(guide=F) +
  scale_x_discrete("Number of heuristics selected",breaks=seq(1,8)) +
  scale_y_continuous("Percentage First Equal") +
  theme_bw(base_size=9)
ggsave(filename="hh_vs_base_top8_fe.pdf", plot=p, width=3, height=3, units="in")
p
p = ggplot(top.N.results, aes(x=n)) + 
  geom_line(aes(y=hhdevmean,linetype="Hyperheuristic")) + 
  geom_point(aes(y=hhdevmean)) +
  geom_line(aes(y=basedevmean,linetype="Baseline")) +
  geom_point(aes(y=basedevmean)) +
  scale_shape_discrete("") +
  scale_linetype_discrete(guide=F) +
  scale_x_discrete("Number of heuristics selected",breaks=seq(1,8)) +
  scale_y_log10("Mean deviation (%)") +
  theme_bw(base_size=9)
ggsave(filename="hh_vs_base_top8_dev.pdf", plot=p, width=3, height=3, units="in")
p
#-------------------------------------------------------------------------------
# If we are splitting out image segmentation, now is the time to analyze it
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
# Here we calculate difficulties and build trees and plots for two heuristics
# as shown in the paper. To produce similar results for other or new heuristics
# simply change the names and adjust the metrics you would like to exclude from
# the analysis.
#-------------------------------------------------------------------------------
# Get difficulty for each instance: deviation / median deviation
results.med.dev <- int.results %>%
  group_by(graphname) %>%
  summarize(meddev=median(dev))
results.diff <- int.results %>%
  left_join(results.med.dev, by="graphname") %>%
  mutate(normdev=dev/meddev) %>%
  mutate(normdev=ifelse(normdev>=10,10,normdev))
#-------------------------------------------------------------------------------
# Reduce metrics set to even more interpretable metrics
more.remove<-c("log_norm_ev1","log_norm_ev2","log_ev_ratio",
               # SKEWS
               "core_skew_positive", "avg_deg_conn_skew_positive",
               "avg_neighbor_deg_skew_positive", "weight_skew_positive",
               "deg_skew_positive", "clust_skew_positive",
               "core_log_abs_skew", "avg_deg_conn_log_abs_skew",
               "avg_neighbor_deg_log_abs_skew", "weight_log_abs_skew",
               "deg_log_abs_skew","clust_log_abs_skew",
               # KURTOSIS
               "core_log_kurtosis","avg_deg_conn_log_kurtosis",
               "avg_neighbor_deg_log_kurtosis","weight_log_kurtosis",
               "deg_log_kurtosis","clust_log_kurtosis")
simple.metrics <- red.metrics[,!(names(red.metrics) %in% more.remove)]
#-------------------------------------------------------------------------------
# Figure 8a: FESTA2002GPR
metrics.removed <-c("log_m","log_n","core_stdev","avg_neighbor_deg_min",
                    "avg_deg_conn_mean","core_mean","avg_deg_conn_min",
                    "avg_neighbor_deg_mean")
simpler.metrics <- simple.metrics[,!(names(simple.metrics) %in% metrics.removed)]
res.sng <- results.diff  %>%
  filter(heuristic == "FESTA2002GPR") %>%
  left_join(simpler.metrics, by="graphname")
t <- rpart(normdev ~ ., data=res.sng[,9:ncol(res.sng)], minbucket=100)
# Tree R^2
1 - sum((predict(t) - res.sng$normdev)^2)/sum((mean(res.sng$normdev) - res.sng$normdev)^2)
# Tree (manually drawn in paper)
prp(t, digits=4, varlen=0, type=1, fallen.leaves=TRUE, tweak=0.8,
    node.fun=function(x, labs, digits, varlen) {
      labs <- lapply(labs, function(l) { substr(l,1,4) })
      return(paste(labs, "\nn=",x$frame$n, sep=""))
    })
# Scatterplot
p <- 
  ggplot(res.sng, aes(x=mis,y=deg_mean,color=normdev)) +
  scale_color_continuous(guide=FALSE, low="blue", high="red") +
  geom_point(size=1) + theme_bw(base_size=9) +
  geom_segment(aes(x=0.1678,y=0.0000, xend=0.1678,yend=1.0000),color="grey",size=0.5) +
  geom_segment(aes(x=0.1678,y=0.0200, xend=1.0000,yend=0.0200),color="grey",size=0.5) +
  annotate("text",x=0.8,y=0.01,label="Easiest (0.2)",color="black",size=2) +
  annotate("text",x=0.8,y=0.05,label="Average (1.2)",color="black",size=2) +
  scale_x_continuous("Max. Indep. Set") +
  scale_y_log10("Density", breaks=c(0.0001,0.001,0.01,0.1,1))
p
ggsave(filename="festa2002_plot.pdf", plot=p, width=3, height=3, units="in")
#-------------------------------------------------------------------------------
# Figure 8b: PALUBECKIS2006
metrics.removed <- c("log_n","log_m","avg_deg_conn_stdev")
simpler.metrics <- simple.metrics[,!(names(simple.metrics) %in% metrics.removed)]
res.sng <- results.diff %>%
  filter(heuristic == "PALUBECKIS2006") %>%
  left_join(simpler.metrics, by="graphname")
t <- rpart(normdev ~ ., data=res.sng[,9:ncol(res.sng)],minbucket=110)
# Tree R^2
1 - sum((predict(t) - res.sng$normdev)^2)/sum((mean(res.sng$normdev) - res.sng$normdev)^2)
# Tree (manually drawn in paper)
prp(t, digits=4, varlen=0, type=1, fallen.leaves=TRUE, tweak=0.8,
    node.fun=function(x, labs, digits, varlen) {
      labs <- lapply(labs, function(l) { substr(l,1,4) })
      return(paste(labs, "\nn=",x$frame$n, sep=""))
    })
# Scatterplot
p <- 
  ggplot(res.sng, aes(x=mis,y=deg_mean,color=normdev)) +
  scale_color_continuous(guide=FALSE, low="blue", high="red") +
  geom_point(size=1) + theme_bw(base_size=9) +
  geom_segment(aes(x=0.3326,y=0.0000, xend=0.3326,yend=1.0000),color="grey",size=0.5) +
  geom_segment(aes(x=0.3326,y=0.0076, xend=1.0000,yend=0.0076),color="grey",size=0.5) +
  geom_segment(aes(x=0.5160,y=0.0000, xend=0.5160,yend=0.0076),color="grey",size=0.5) +
  annotate("text",x=0.5,y=0.05,label="Hardest (4.5)",color="black",size=2) +
  annotate("text",x=0.8,y=0.0008,label="Harder (3.5)",color="black",size=2) +
  annotate("text",x=0.43,y=0.00013,label="Hard (1.6)",color="black",size=2) +
  scale_x_continuous("Max. Indep. Set") +
  scale_y_log10("Density", breaks=c(0.0001,0.001,0.01,0.1,1))
p
ggsave(filename="pal2006_plot.pdf", plot=p, width=3, height=3, units="in")
#-------------------------------------------------------------------------------


################################################################################
# Section 5.2: Understanding heuristic ideas
# Here we provide a a generic function to create the density plots as shown
# in the paper. By changing the heuristic features data, or by adding new
# heuristics/instances, the results can be updated
#-------------------------------------------------------------------------------
# See what labels look like they might make sense
res.class <- results %>%
  filter(rank2 == 1.0,
         graphname %in% int.graphs$graphname) %>%
  left_join(heuristics, by="heuristic")
# The desireable property is that one feature is not dominant over the others
table(res.class$classification)  # Iterated Local Search
table(res.class$Hybrid)  # Balanced
table(res.class$Initialization)  # Balanced
table(res.class$Population)  # Balanced
table(res.class$Evolutionary)  # No vs 1+2Parent
table(res.class$Memory)  # Balanced
table(res.class$Perturbation)  # Balanced
table(res.class$problem)  # Not very balanced
#-------------------------------------------------------------------------------
# Define the plotting function
nnplot <- function(xdata, ydata, outcome, pts, filename, k, rad, xlab, ylab, title) {
  df <- data.frame(x=xdata, y=ydata, outcomes=outcome)
  
  # Determine area
  xmin <- ifelse(min(df$x) < 0, min(df$x) * 1.01, min(df$x) * 0.99)
  xmax <- ifelse(max(df$x) < 0, max(df$x) * 0.99, max(df$x) * 1.01)
  ymin <- ifelse(min(df$y) < 0, min(df$y) * 1.01, min(df$y) * 0.99)
  ymax <- ifelse(max(df$y) < 0, max(df$y) * 0.99, max(df$y) * 1.01)
  sds <- c((xmax - xmin)^2, (ymax - ymin)^2)
  
  # Setup grid
  xvals <- seq(xmin, xmax, length.out=pts)
  yvals <- seq(ymin, ymax, length.out=pts)
  points <- as.matrix(expand.grid(xvals, yvals))
  nn <- seq(0, 0, length.out=pts*pts)
  min.dist <- seq(100, 100, length.out=pts*pts)
  n.range <- seq(0, 0, length.out=pts*pts)
  
  # Calculate nearest neighbors
  for (i in 1:nrow(points)) {
    px <- points[i,1]
    py <- points[i,2]
    df.dist <- df %>%
      mutate(dist=sqrt((x-px)^2+(y-py)^2)) %>%
      arrange(dist)
    nn[i] <- sum(df.dist[1:k,3] == 1) / k
    n.range[i] <- sum(df.dist[1:k,4] < rad)
  }
 
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
#   axis(1, at=c(1,2,3,4,5), labels=c(expression(10^1),
#                                     expression(10^2),
#                                     expression(10^3),
#                                     expression(10^4),
#                                     expression(10^5)))
#   axis(2, at=c(2,3,4,5,6,7), labels=c(expression(10^2),
#                                   expression(10^3),
#                                   expression(10^4),
#                                   expression(10^5),
#                                   expression(10^6),
#                                   expression(10^7)))
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
  left_join(simple.metrics, by="graphname") %>%
  filter(weight_mean>=-3.5e-01)
nnplot(log10(exp(res.EVO$log_n)), log10(exp(res.EVO$log_m)), as.numeric(res.EVO$Evolutionary=="Yes"),
       400, "nn_evo.png", k=20, rad=1, xlab="Number of Nodes", ylab="Number of Edges", title="Evolutionary Algorithm")
#-------------------------------------------------------------------------------
# MEMORY (TABU SEARCH)
res.MEM <- res.class %>%
  select(graphname,Memory) %>%
  mutate(Memory=as.factor(Memory)) %>%
  left_join(simple.metrics, by="graphname") %>%
  filter(weight_mean>=-3.5e-01)
nnplot(log10(exp(res.MEM$log_n)), log10(exp(res.MEM$log_m)), as.numeric(res.MEM$Memory=="Yes"),
       400, "nn_mem.png", k=20, rad=1, xlab="Number of Nodes", ylab="Number of Edges", title="Tabu Search")
#-------------------------------------------------------------------------------