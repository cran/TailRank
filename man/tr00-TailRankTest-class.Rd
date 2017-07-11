\name{TailRankTest-class}
\docType{class}
\alias{TailRankTest-class}
\alias{getStatistic}

\title{The TailRankTest Class}
\description{This is the class representation for the results of a
  tail-rank test to find biomarkers in a microarray data set. It
  includes methods for summarizing and plotting the results of the
  test.}
\section{Creating objects}{
  Although objects can be created, as usual, using \code{new}, the only
  reliable way to create valid objects is to use the
  \code{\link{TailRankTest}} function. See the description of that
  function for details on how the tail-rank test works.
}
\section{Slots}{
  \describe{
    \item{\code{statistic}:}{a numeric vector containng the tail-rank
      statistic for each row (gene) in a microarray data set}
    \item{\code{direction}:}{a character string representing the
      direction of the test; can be "up", "down", or "two-sided"}
    \item{\code{N1}:}{an integer; the numnber of samples in the "base" or
      "healthy" group}
    \item{\code{N2}:}{an integer; the number of samples in the "test" or
      "cancer" group}
    \item{\code{specificity}:}{a real number between 0 and 1; the desired
      specificity used in the test to estimate a quantile from the
      "base" group}
    \item{\code{tolerance}:}{a real number between 0 and 1; the upper
      tolerance bound used to estimate the threshold}
    \item{\code{confidence}:}{a real number between 0 and 1; the
      confidence level that there are no false positives}
    \item{\code{cutoff}:}{an integer; the maximum expected value of the
      statistic under the null hypothesis}
    \item{\code{model}:}{a character string describing the model
      (binomial or beta-binomial) used to decide on cutoffs for significance}
    \item{\code{tau}:}{a numeric vector or NULL; gene-by-gene upper
      bounds for significance}
    \item{\code{rho}:}{a numeric vector or NULL; gene-by-gene lower
      bounds for significance}
  }
}
\section{Methods}{
  \describe{
    \item{summary(object, ...)}{Display a summary of the TailRankTest \code{object}}
    \item{hist(x, overlay, ...)}{Plot a histogram of the statistic in
      the TailRankTest object \code{x}. The optional argument
      \code{overlay} is a logical flag. If \code{overlay=TRUE}, then the
      histogram is overlain with a curve representing the null
      distribution. The default value of \code{overlay} is \code{FALSE}.}
    \item{as.logical(x, ...)}{Convert the TailRankTest object
      \code{x} into a logical vector, which takes on a \code{TRUE} value
      whenever the tail-rank statistic exceeds the significance cutoff.}
    \item{getStatistic(object, ...)}{Obtain the vector of tail-rank
      statistics contained in \code{object}.}
  }
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\seealso{\code{\link{TailRankTest}},
  \code{\link{tailRankPower}},
  \code{\link{biomarkerPowerTable}},
  \code{\link[oompaBase]{matrixMean}},
  \code{\link{toleranceBound}}
}
\examples{
# generate some fake data to use in the example
nr <- 40000
nc <- 110
fake.data <- matrix(rnorm(nr*nc), ncol=nc)
fake.class <- rep(c(TRUE, FALSE), c(40, 70))

# perform the tail-rank test
null.tr <- TailRankTest(fake.data, fake.class)

# get a summary of the results
summary(null.tr)

# plot a histogram of the statistics
hist(null.tr, overlay=TRUE)

# get the actual statistics
stats <- getStatistic(null.tr)

# get a vector that selects the "positive" calls for the test
is.marker <- as.logical(null.tr)

# the following line should evaluate to the number of rows, nr = 40000
sum( is.marker == (stats > null.tr@cutoff) )
}
\keyword{univar}
\keyword{multivariate}
\keyword{htest}
\keyword{nonparametric}
