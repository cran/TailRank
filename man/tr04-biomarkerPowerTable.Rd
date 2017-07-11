\name{biomarkerPowerTable}
\alias{biomarkerPowerTable}
\title{Power tables for the tail-rank test}
\description{
  Compute an array of power tables for the tail-rank.test.
}
\usage{
biomarkerPowerTable(G, N1=20, N2=seq(25, 250, by=25),
                    psi = c(0.95, 0.99), conf=0.99,
                    phi = seq(0.10, 0.50, by = 0.05), model="bb")
}
\arguments{
  \item{G}{An integer; the number of genes being assessed as potential
    biomarkers. Statistically, the number of hypotheses being tested.}
  \item{N1}{An integer; the number of "train" or "healthy" samples used.}
  \item{N2}{An integer; the number of "test" or "cancer" samples used.}
  \item{psi}{A real number between 0 and 1; the desired specificity of
    the test.}
  \item{conf}{A real number between 0 and 1; the confidence level of the
    results. Can be obtained by subtracting the family-wise Type I error
    from 1.}
  \item{phi}{A real number between 0 and 1; the sensitivity that one
    would like to be able to detect, conditional on the specificity.}
  \item{model}{A character string that determines whether power and
    significance are computed from abinomial or a beta-binomial (bb)
    model.}
}
\value{
  Returns a list of objects of the \code{\link{BMPT}} class. Each item
  in the list consists of a two-dimensional table (indexed by the sample
  sizes \code{N} and the sensitivities \code{phi}) with scalars
  recording the values of \code{G}, \code{conf}, and \code{psi} that
  were used to generate it.
}

\note{Default values of the optional arguments (\code{N}, \code{psi},
  \code{conf}, \code{phi})are included in the usage examples.
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\seealso{\code{\link{TailRankTest}},
  \code{\link{tailRankPower}},
  \code{\link{biomarkerPowerTable}},
  \code{\link[oompaBase]{matrixMean}},
  \code{\link{toleranceBound}}
}
\examples{
stuff <- biomarkerPowerTable(10000, 20,
                               c(10, 20, 50, 100, 250, 500),
                               c(0.95, 0.99),
                               c(0.99, 0.95),
                               seq(0.1, 0.7, by=0.1))
lapply(stuff, summary)
}
\keyword{univar}
\keyword{multivariate}
\keyword{htest}
\keyword{nonparametric}
