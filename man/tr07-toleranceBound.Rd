\name{toleranceBound}
\alias{toleranceBound}
\title{Upper tolerance bounds on normal quantiles}
\description{
The function \code{toleranceBound} computes theoretical upper tolerance
bounds on the quantiles of the standard normal distribution. These can
be used to produce reliable data-driven estimates of the quantiles in
any normal distribution.
}
\usage{
toleranceBound(psi, gamma, N)
}
\arguments{
  \item{psi}{A real number between 0 and 1 giving the desired quantile }
  \item{gamma}{A real number between 0 and 1 giving the desired
    tolerance bound}
  \item{N}{An integer giving the number of observations used to estimate
  the quantile}
}
\details{

Suppose that we collect \eqn{N} observations from a normal distribution
with unknown mean and variance, and wish to estimate the \eqn{95}th
percentile of the distribution. A simple point estimate is given by
\eqn{\tau = \bar{X} + 1.68s}. However, only the mean of the distribution is
less than this value \eqn{95\%} of the time.  When \eqn{N=40}, for example,
almost half of the time (\eqn{43.5\%}), fewer than \eqn{95\%} of the
observed values will be less than \eqn{\tau}. This problem is addressed by
constructing a statistical tolerance interval (more precisely, a one-sided
tolerance bound) that contains a given fraction, \eqn{\psi}, of the
population with a given confidence level, \eqn{\gamma} [Hahn and Meeker,
1991]. With enough samples, one can obtain distribution-free tolerance
bounds [op.\ cit., Chapter 5]. For instance, one can use bootstrap or
jackknife methods to estimate these bounds empirically.

Here, however, we assume that the measurements are normally distributed. We
let \eqn{\bar{X}} denote the sample mean and let \eqn{s} denote the sample
standard deviation. The upper tolerance bound that, \eqn{100 \gamma\%} of
the time, exceeds \eqn{100 \psi\%} of \eqn{G} values from a normal
distribution is approximated by \eqn{X_U = \bar{X} + k_{\gamma,\psi}s},
where
\deqn{
  k_{\gamma, \psi} = {z_{\psi} + \sqrt{z_{\psi}^2 - ab} \over a},
}\deqn{
  a = 1-{z_{1-\gamma}^2\over 2N-2},
}\deqn{
  b = z_{\psi}^2 - {z_{1-\gamma}^2\over N},
}
and, for any \eqn{\pi}, \eqn{z_\pi} is the critical value of the normal
distribution that is exceeded with probability \eqn{\pi} [Natrella, 1963].
}
\value{
  Returns the value of \eqn{k_{\gamma, \psi}} with the property that the
  \eqn{\psi}th quantile will be less than the estimate \eqn{X_U =
    \bar{X} + k_{\gamma,\psi}s} (based on \eqn{N} data points) at least
  \eqn{100 \gamma\%} of the time. 
}
\references{
Natrella, M.G. (1963) \emph{Experimental Statistics}.
NBS Handbook 91, National Bureau of Standards, Washington DC.

Hahn, G.J. and Meeker, W.Q. (1991)
\emph{Statistical Intervals: A Guide for Practitioners}.
John Wiley and Sons, Inc., New York.
}
\note{Lower tolerance bounds on quantiles with \code{psi} less than
  one-half can be obtained as \eqn{X_U = \bar{X} - k_{\gamma,1-\psi}s},}
\author{Kevin R. Coombes <krc@silicovore.com>}
\examples{
N <- 50
x <- rnorm(N)
tolerance <- 0.90
quant <- 0.95
tolerance.factor <- toleranceBound(quant, tolerance, N)

# upper 90% tolerance bound for 95th percentile
tau <- mean(x) + sd(x)*tolerance.factor

# lower 90% tolerance bound for 5th percentile
rho <- mean(x) - sd(x)*tolerance.factor

# behavior of the tolerance bound as N increases
nn <- 10:100
plot(nn, toleranceBound(quant, tolerance, nn))

# behavior of the bound as the tolerance varies
xx <- seq(0.5, 0.99, by=0.01)
plot(xx, toleranceBound(quant, xx, N))
}
\keyword{univar}
\keyword{htest}
