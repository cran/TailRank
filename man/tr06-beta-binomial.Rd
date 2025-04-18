\name{BetaBinomial}
\alias{dbb}
\alias{pbb}
\alias{qbb}
\alias{rbb}
\title{The Beta-Binomial Distribution}
\description{
Density, distribution function, quantile function, and random generation
for the beta-binomial distribution.  A variable with a beta-binomial
distribution is distributed as binomial distribution with parameters
\code{N} and \code{p}, where the probability \code{p} of success iteself
has a beta distribution with parameters \code{u} and \code{v}.
}
\usage{
dbb(x, N, u, v, log = FALSE)
pbb(q, N, u, v)
qbb(p, N, u, v)
rbb(n, N, u, v)
}
\arguments{
  \item{x}{vector of qauntiles}
  \item{q}{vector of quantiles}
  \item{p}{vector of probabilities}
  \item{n}{number of observations}
  \item{N}{number of trials ( a positive integer)}
  \item{u}{first positive parameter of the beta distribution}
  \item{v}{second positive parameter of the beta distribution}
  \item{log}{A logical value; if true, values are returned on the log scale}
}
\details{
  The beta-binomial distribution with parameters \eqn{N}, \eqn{u}, and
  \eqn{v} has density given by
  \deqn{
    choose(N, x) * Beta(x + u, N - x + v) / Beta(u,v)
  }
  for \eqn{u > 0}, \eqn{v > 0}, a positive integer \eqn{N}, and any
  nonnegative integer \eqn{x}. Although one can express the integral
  in closed form using generalized hypergeometric functions, the
  implementation of  distribution function used here simply relies on the
  the cumulative sum of the density.

  The mean and variance of the beta-binomial distribution can be
  computed explicitly as
  \deqn{
    \mu = \frac{Nu}{u+v}
  }
  and
  \deqn{
    \sigma^2 = \frac{Nuv(N+u+v)}{(u+v)^2 (1+u+v)}
  }
}
\value{
  \code{dbb} gives the density, \code{pbb} gives the distribution function,
  \code{qbb} gives the quantile function, and \code{rbb} generates random
  deviates. 
}
\seealso{ \code{\link[stats]{dbeta}} for the beta distribution and
  \code{\link[stats]{dbinom}} for the binomial distribution.
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\examples{
# set up parameters
w <- 10
u <- 0.3*w
v <- 0.7*w
N <- 12
# generate random values from the beta-binomial
x <- rbb(1000, N, u, v)

# check that the empirical summary matches the theoretical one
summary(x)
qbb(c(0.25, 0.50, 0.75), N, u, v)

# check that the empirpical histogram matches te theoretical density
hist(x, breaks=seq(-0.5, N + 0.55), prob=TRUE)
lines(0:N, dbb(0:N, N, u,v), type='b')
}
\keyword{univar}
\keyword{htest}
