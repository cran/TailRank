
# TR-tables.R
# Copyright, Kevin R. Coombes, 2004.

# This file both gives an example of how to use the tail-rank-power
# functions and genrates the actual table sused in the paper on the
# tail-rank statistic.

#######################
# First we generate Table 1 of the paper. This table contains the
# maximum expected value of the statistic over G independent tests.

# Left half of table 1, with univariate specificity psi = 0.99
psi.0 <- 0.99
confide <- rev(c(0.8, 0.95, 0.99))
ng <- c(100, 1000, 10000, 100000)
ns <- c(10, 20, 50, 100, 250, 500)

formal.cut <- array(0, c(length(ns), length(ng), length(confide)))
for (i in 1:length(ng)) {
  for (j in 1:length(ns)) {
    formal.cut[j, i, ] <- tail.rank.cutoff(ng[i], ns[j], psi.0, confide)
  }
}
dimnames(formal.cut) <- list(ns, ng, confide)
formal.cut

# Right half of table 1, with univariate specificity psi = 0.95
psi.0 <- 0.95
formal.cut <- array(0, c(length(ns), length(ng), length(confide)))
for (i in 1:length(ng)) {
  for (j in 1:length(ns)) {
    formal.cut[j, i, ] <- tail.rank.cutoff(ng[i], ns[j], psi.0, confide)
  }
}
dimnames(formal.cut) <- list(ns, ng, confide)
formal.cut

#######################
# Now we get ready to make the second table
stuff <- biomarker.power.table(10000,
                               c(10, 20, 50, 100, 250, 500),
                               c(0.95, 0.99),
                               c(0.99, 0.95),
                               seq(0.1, 0.7, by=0.1))

# Here are the four panels of table 2
lapply(stuff, summary)


## cleanup
rm(psi.0, confide, ng, ns, formal.cut, stuff)
