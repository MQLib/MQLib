library(dplyr)
library(ggplot2)

# Load graph sizes
graph.sizes <- read.csv("metrics.csv",stringsAsFactors=F) %>%
  mutate(n = exp(log_n)) %>%
  select(filename, n)

# Load data for a full 20 minute run of 60 randomly selected instances with
# all 37 heuristics
rt.data <- read.csv("runtime_data.csv", stringsAsFactors=F) %>%
  filter(runtime <= 1200) %>%
  left_join(graph.sizes, by=c("graphname"="filename"))

# Add a gap column that, for each instance-heuristic pair, shows the gap from
# the best solution FOR THAT HEURISTIC for each of the intermediate solutions
rt.data <- rt.data %>%
  group_by(graphname, heuristic) %>%
  mutate(gap = 1 - (objective / max(objective))) %>%
  ungroup()

# Summarize the data to find the last arrival time for each pair of
# heuristic and instance
last.arrivals <- rt.data %>%
  group_by(graphname, heuristic) %>%
  summarize(runtime = max(runtime)) %>%
  left_join(graph.sizes, by=c("graphname"="filename"))

# Compute coverage and runtime over a large grid of n coeffs
n.coefs <- seq(0.2, 1.0, 0.05)
min.val <- 120   # The minimum runtime limit
max.val <- 1200  # The maximum

# Return total runtime and coverage
run.info <- do.call(rbind, lapply(n.coefs, function(k) {
  print(k)
  # Calculate the runtimes this coefficient would lead to for the
  # 60 instance subset:
  rt.60s <- pmin(pmax(k * last.arrivals$n, min.val), max.val)
  # All instances
  rt.all <- pmin(pmax(k * graph.sizes$n, min.val), max.val)
  # Determine the gaps this would have resulted in (on the subset)
  gaps <- rt.data %>%
    mutate(rtlimit = pmin(pmax(k * n, min.val), max.val)) %>%
    group_by(graphname, heuristic) %>%
    filter(runtime <= rtlimit) %>%
    # No solution is the same as solution 0, with gap 1
    summarize(mingap = ifelse(n() > 0, min(gap), 1))
  return(data.frame(k=k,
             days=sum(rt.all)*37/3600/24/60,
             cost=sum(rt.all)*37/3600*.07,
             coverage=mean(rt.60s >= last.arrivals$runtime),
             gap=mean(gaps$mingap)))
}))
run.info

# Plot various metrics versus the coefficient
ggplot(run.info, aes(x=k, y=coverage)) + geom_line()
ggplot(run.info, aes(x=k, y=gap)) + geom_line()
ggplot(run.info, aes(x=k, y=cost)) + geom_line()