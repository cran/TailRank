\name{tailRankPower}
\alias{tailRankPower}
\alias{tailRankCutoff}
\title{Power of the tail-rank test}
\description{
  Compute the significance level and the power of a tail-rank test.
}
\usage{
tailRankPower(G, N1, N2, psi, phi, conf = 0.95,
              model=c("bb", "betabinom", "binomial"))
tailRankCutoff(G, N1, N2, psi, conf,
               model=c("bb", "betabinom", "binomial"),
               method=c('approx', 'exact'))
                           
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{G}{An integer; the number of genes being assessed as potnetial
    biomarkers. Statistically, the number of hypotheses being tested.}
  \item{N1}{An integer; the number of "train" or "healthy" samples used.}
  \item{N2}{An integer; the number of "test" or "cancer" samples used.}
  \item{psi}{A real number between 0 and 1; the desired specificity of
    the test.}
  \item{phi}{A real number between 0 and 1; the sensitivity that one
    would like to be able to detect, conditional on the specificity.}
  \item{conf}{A real number between 0 and 1; the confidence level of the
    results. Can be obtained by subtracting the family-wise Type I error
    from 1.}
  \item{model}{A character string that determines whether significance
    and power are computed based on a binomial or a beta-binomial (bb)
    model.}
  \item{method}{A character string; either "exact" or "approx". The
    deafult is to use a Bonferroni approximation.}
}
\details{
A power estimate for the tail-rank test can be obtained as follows.
First, let X ~ Binom(N,p) denote a binomial random variable. Under the
null hypotheis that cancer is not different from normal, we let
\eqn{p = 1 - \psi} be the expected  proportion of successes in a test of
whether the value exceeds the psi-th quantile.  Now let
\deqn{\alpha = P(X > x,| N, p)}
be one such binomial measurement.  When we make \eqn{G} independent
binomial measurements, we take
\deqn{conf = P(all\ G\ of\ the\ X's \le x | N, p).}
(In our paper on the tail-rank statistic, we write everything in terms
of \eqn{\gamma = 1 - conf}.) Then we have
\deqn{conf = P(X \le x | N, p)^G = (1 - alpha)^G.}
Using a Bonferroni-like approximation, we can take
\deqn{conf ~= 1 - \alpha*G.}
 Solving for \eqn{\alpha}, we find that
\deqn{\alpha ~= (1-conf)/G.}
So, the cutoff that ensures that in multiple experiments, each looking
at \eqn{G} genes in \eqn{N} samples, we have confidence level \eqn{conf}
(or significance level \eqn{\gamma = 1 - conf}) of no false positives is
computed by the function \code{tailRankCutoff}.

The final point to note is that the quantiles are also defined
in terms of \eqn{q = 1 - \alpha}, so there are lots of disfiguring "1's"
in the implementation.

Now we set \eqn{M} to be the significance cutoff using the procedure
detailed above.  A gene with sensitivity \eqn{\phi} gets detected if the
observed number of cases above the threshold is greater than or equal to
\eqn{M}. The \code{tailRankPower} function implements formula (1.3) of
our paper on the tail-rank test.

}
\value{
  \code{tailRankCutoff} returns an integer that is the
  maximum expected value of the tail rank statistic under the null
  hypothesis.
    
  \code{tailRankPower} returns a real numbe between 0 and 1
  that is the power of the tail-rank test to detect a marker with
  true sensitivity equal to \eqn{phi}.
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\seealso{\code{\link{TailRankTest}},
  \code{\link{tailRankPower}},
  \code{\link{biomarkerPowerTable}},
  \code{\link[oompaBase]{matrixMean}},
  \code{\link{toleranceBound}}
}
\examples{
psi.0 <- 0.99
confide <- rev(c(0.8, 0.95, 0.99))
nh <- 20
ng <- c(100, 1000, 10000, 100000)
ns <- c(10, 20, 50, 100, 250, 500)
formal.cut <- array(0, c(length(ns), length(ng), length(confide)))
for (i in 1:length(ng)) {
  for (j in 1:length(ns)) {
    formal.cut[j, i, ] <- tailRankCutoff(ng[i], nh, ns[j], psi.0, confide)
  }
}
dimnames(formal.cut) <- list(ns, ng, confide)
formal.cut

phi <- seq(0.1, 0.7, by=0.1)
N <- c(10, 20, 50, 100, 250, 500)
pows <- matrix(0, ncol=length(phi), nrow=length(N))
for (ph in 1:length(phi)) {
  pows[, ph] <-  tailRankPower(10000, nh, N, 0.95, phi[ph], 0.9)
}
pows <- data.frame(pows)
dimnames(pows) <- list(as.character(N), as.character(round(100*phi)))
pows
}
\keyword{univar}
\keyword{multivariate}
\keyword{htest}
\keyword{nonparametric}
