# Copyright (C) Kevin R. Coombes, 2007-2012

# tail-rank-aux.R
# Copyright, Kevin R. Coombes, 2004

############################
# Finally, we include a function to include tolerance bounds on the
# quantiles of a normal distribution. The reference for this function
# is the National Bureau of Standards handbook by Natrella, 1963.
toleranceBound <- function(psi, gamma, N) {
  zg <- qnorm(1-gamma)
  zp <- qnorm(psi)
  a <- 1 - zg^2/(2*N-2)
  b <- zp^2 - zg^2/N
  k1 <- (zp + sqrt(zp^2-a*b))/a
  k1
}
