# Copyright (C) Kevin R. Coombes, 2007-2012

dbb <- function(x, N, u, v, log = FALSE) {
  logval <- lbeta(x+u, N-x+v) - lbeta(u,v) + lchoose(N,x)
  if (log) {
    ret <- logval
  } else {
    ret <- exp(logval)
  }
  ret
}

pbb <- function(q, N, u, v) {
  sapply(q, function(xx) sum(dbb(0:xx, N, u, v)))
}

qbb <- function(p, N, u, v) {
  pp <- cumsum(dbb(0:N, N, u, v))
  sapply(p, function(x) sum(pp < x))
}

rbb <- function(n, N, u, v) {
  p <- rbeta(n, u, v)
  rbinom(n, N, p)
}
