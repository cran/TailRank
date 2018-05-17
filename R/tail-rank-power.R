# Copyright (C) Kevin R. Coombes, 2007-2012

# tail-rank-power.R
# Copyright, Kevin R. Coombes, 2004.

require('methods')

#########################################################################
# A power estimate can be obtained as follows.  First, let X ~ Binom(N,p)
# denote a binomial random variable. Under the null hypotheis that cancer
# is not different from normal, we can choose p = 1 - psi for the expected
# proportion of successes in a test of whether the value exceeds the
# psi-th quantile.  Now let
#	alpha = P(X > x,| N, p)
# for one such binomial measurement.  When we make G independent binomial
# measurements, we take
#	conf = P(all G of the X's <= x | N, p).

# (In our paper on the tail-rank statistic, we write everything in terms
# of gamma = 1 - conf.) Then we have
#	conf = P(X <= x | N, p)^G = (1 - alpha)^G.
# Using a Bonferroni-like approximation, we can take
#	conf ~= 1 - # alpha*G.
# Solving for alpha, we find that
#	alpha ~= (1-conf)/G.
# So, the cutoff that ensures that in multiple experiments, each looking
# at G genes in N samples, we have confidence level conf (or significance
# level gamma = 1 - conf) of no false positives is given by the following
# function. The final point to note is that the quantiles are also defined
# in terms of q = 1 - alpha, so there are lots of disfiguring "1's" in the
# implementation.
tailRankCutoff <- function(G, N1, N2, psi, conf, 
                           model=c("bb", "betabinom", "binomial"),
                           method=c('approx', 'exact')) {
  method <- match.arg(method)
  if (method == "approx") {
    q <- 1 - (1-conf)/G
  } else {
    q <- conf^(1/G)
  }
  model <- match.arg(model)
  if(model =="binomial") {
    value <- qbinom(q, N2, 1-psi)
  } else {
    W <- N1+2
    value <- sapply(N2, function(N) qbb(q, N, W*(1-psi), W*psi))
  }
  value
}

# Now we again set M to be the significance cutoff using the procedure
# detailed above.  A gene with sensitivity phi gets detected if the
# observed number of cases above the threshold is greater than or equal
# to M. So, we have a function that implements formula (1.3) of the paper
# on the tail-rank test.
tailRankPower <- function(G, N1, N2, psi, phi, conf=0.95,
                          model=c("bb", "betabinom", "binomial")) {
  model <- match.arg(model)
  M <- tailRankCutoff(G, N1, N2, psi, conf, model=model)
  if(model =="binomial") {
    value <- 1 - pbinom(M, N2, phi)
  } else {
    W <- N1+2
    value <- sapply(1:length(N2), function(i) 1 - pbb(M[i], N2[i], W*phi, W*(1 - phi)))
  }
  value
}

#########################################################################
# Now we introduce a class to generate power tables.

setClass('BMPT',
         slots = c(G='numeric',
                   psi='numeric',
                   conf='numeric',
                   power='data.frame')
         )


BMPT <- function(G, psi, conf, power) {
  new('BMPT', G=G, psi=psi, conf=conf, power=power)
}

setMethod('summary', signature(object='BMPT'),
          function(object, ...) {
              object
          })

setMethod('print', signature(x='BMPT'),
          function(x, ...) {
  cat(paste('G = ', x@G, ', psi = ', x@psi,
            ', conf = ', x@conf, '\nPower:\n', sep=''))
  print(x@power)
})

biomarkerPowerTable <-
  function(G,
           N1=20,
           N2=seq(25, 250, by=25),
           psi=c(0.95, 0.99),
           conf=0.99,
           phi=seq(0.10, 0.50, by=0.05),
           model="bb") {
  val <- list()
  i <- 1
  for (g in G) {
    for (p in psi) {
      for (gam in conf) {
	pows <- matrix(0, ncol=length(phi), nrow=length(N2))
	for (ph in 1:length(phi)) {
          pows[, ph] <-  tailRankPower(g, N1, N2, p, phi[ph], gam, model)
	}
	pows <- data.frame(pows)
	dimnames(pows) <- list(as.character(N2), as.character(round(100*phi)))
        extra <- BMPT(g, p, gam, pows)
        val[[i]] <-  extra
        i <- i + 1
      }
    }
  }
  val
}

