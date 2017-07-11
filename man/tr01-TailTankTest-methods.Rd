\name{TailRankTest-methods}
\docType{methods}
\alias{getStatistic,TailRankTest-method}
\alias{summary,TailRankTest-method}
\alias{hist,TailRankTest-method}
\alias{as.logical,TailRankTest-method}

\title{Methods for TailRankTest objects}
\description{
  This file describes the methods for an object of the class
  \code{\link[=TailRankTest-class]{TailRankTest}} class.
}
\usage{
\S4method{summary}{TailRankTest}(object, ...)
\S4method{hist}{TailRankTest}(x, overlay = FALSE,
  xlab = "tail-rank statistic", main = "", ...)
\S4method{as.logical}{TailRankTest}(x, ...)
\S4method{getStatistic}{TailRankTest}(object,...)
}
\arguments{
  \item{x}{A \code{\link[=TailRankTest-class]{TailRankTest} object}}
  \item{object}{A \code{\link[=TailRankTest-class]{TailRankTest} object}}
  \item{overlay}{An optional logical flag; defaults to \code{FALSE}.}
  \item{xlab}{A character string}
  \item{main}{A character string}
  \item{\dots}{Extra graphical parameters}
}
\value{
  \item{\code{as.logical}}{Returns a logical vector. \code{TRUE} values
    pick out candidate biomarkers where the tail-rank test statistic
    exceeds the significance cutoff.}
  \item{\code{getStatistic}}{Returns the vector of tail-rank statistics
    contained in \code{object}.}
  \item{\code{hist}}{Invisibly returns the TailRankTest object.}
  \item{\code{summary}}{Invisibly returns the TailRankTest object.}
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\seealso{\code{\link{TailRankTest-class}},
  \code{\link{TailRankTest}},
  \code{\link{tailRankPower}}
}
\examples{
# generate some fake data to use in the example
nr <- 40000
nc <- 110
fake.data <- matrix(rnorm(nr*nc), ncol=nc)
fake.class <- rep(c(TRUE, FALSE), c(40, 70))

# build an object
null.tr <-  TailRankTest(fake.data, fake.class)

# summarize the object
summary(null.tr)

# plot a histogram
hist(null.tr)
hist(null.tr, breaks=70, col='blue', overlay=TRUE)

# get a logical vector that can select those markers
# identified by the test
selector <- as.logical(null.tr)
}
\keyword{manip}
