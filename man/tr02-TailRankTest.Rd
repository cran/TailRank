\name{TailRankTest}
\alias{TailRankTest}

\title{The Tail-Rank Test}
\description{
  Perform a tail-rank test to find candidate biomarkers in a microarray
  data set.
}
\usage{
TailRankTest(data, classes, specificity = 0.95, tolerance = 0.50,
             model=c("bb", "betabinomial", "binomial"),
             confidence = 0.95, direction = "up")
}
\arguments{
  \item{data}{A matrix or data.frame containing numerical measurements on
    which to perform the tail-rank test.}
  \item{classes}{A logical vector or factor splitting the data into two
    parts. The length of this vector should equal the number of columns
    in the \code{data}. The \code{TRUE} portion (or the first level of
    the factor) represents a "base" or "healthy" group of samples; the
    other samples are the "test" or "cancer" group.} 
  \item{specificity}{a real number between 0 and 1; the desired
    specificity used in the test to estimate a quantile from the
    "base" group. This is an optional argument with default value 0.95.}
  \item{tolerance}{a real number between 0 and 1; the upper tolerance
    bound used to estimate the threshold. This is an optional argument
    with default value 0.90.}
  \item{model}{a character string that determines whther significance
    comes from a binomial model or a beta-binomial (bb) model.}
  \item{confidence}{a real number between 0 and 1; the confidence level
    that there are no false positives. This is an optional argument with
    default value 0.50, which is equivalent to ignoring the tolerance.} 
  \item{direction}{a character string representing the direction of the
    test; can be "up", "down", or "two-sided". The default value is "up".} 
}
\details{
This function computes the tail rank statistic for each gene (viewed as
one row of the data matrix). The data is split into two groups. The
first ("base") group is used to estimate a tolerance bound (defaults to
50\%) on a specific quantile (defaults to 95\%) of the distribution of
each gene. The tail-rank statistic is the defined as the number of
samples in the second ("test") group that lie outside the bound. The
test can be applied in the "up", "down", or "two-sided" direction,
depending on the kinds of markers being sought. Also computes the cutoff
for significance based on a confidence level that is "1 - FWER" for a
desired family-wise error rate.
}
\value{
  The return value is an object of class
  \link[=TailRankTest-class]{TailRankTest}.
}
\references{http://bioinformatics.mdanderson.org}
\author{Kevin R. Coombes <krc@silicovore.com>}
\seealso{\code{\link{TailRankTest-class}},
  \code{\link{tailRankPower}},
  \code{\link{biomarkerPowerTable}},
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
