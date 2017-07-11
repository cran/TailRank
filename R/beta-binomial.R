# Copyright (C) Kevin R. Coombes, 2007-2012

dbb <- function(x, N, u, v) {
  beta(x+u, N-x+v)/beta(u,v)*choose(N,x)
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
