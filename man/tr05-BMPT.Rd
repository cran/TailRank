\name{BMPT-class}
\docType{class}
\alias{BMPT-class}
\alias{summary,BMPT-method}
\alias{print,BMPT-method}
\alias{BMPT}

\title{The BMPT Class}
\description{
  A class for producing BioMarker Power Tables (BMPT), and methods for
  accessing them. This class is primarily an implementation detail for
  the function \code{\link{biomarkerPowerTable}}.
}
\section{Creating objects}{
  Although objects can be created using \code{new}, the preferred method
  is to use the constructor function \code{BMPT}. In practice, these
  objects are most likely to be created using the more general interface
  through \code{\link{biomarkerPowerTable}}.
}
\usage{
BMPT(G, psi, conf, power)
\S4method{print}{BMPT}(x,...)
\S4method{summary}{BMPT}(object,...)
}
\arguments{
  \item{G}{A positive integer.}
  \item{psi}{A real number between 0 and 1.}
  \item{conf}{A real number between 0 and 1.}
  \item{power}{A data frame.}
  \item{x}{A \code{BMPT} object.}
  \item{object}{A \code{BMPT} object.}
  \item{\dots}{Extra graphical parameters}
}
\section{Slots}{
  \describe{
  \item{\code{G}:}{A positive integer; the number of genes being
    assessed as potential biomarkers. Statistically, the number of
    hypotheses being tested.}
  \item{\code{psi}:}{A real number between 0 and 1; the desired
    specificity of the test.}
  \item{\code{conf}:}{A real number between 0 and 1; the confidence
    level of the results. Can be obtained by subtracting the family-wise
    Type I error from 1.}
  \item{\code{power}:}{A data frame containing the power
    computations. The rows are indexed by the sample size and the
    columns by the sensitivity.}
  }
}
\section{Methods}{
  \describe{
    \item{print(x, ...)}{Print the power table \code{x}.}
    \item{summary(object, ...)}{Summarize the power table \code{object}.}
  }
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\seealso{\code{\link{TailRankTest}},
  \code{\link{tailRankPower}},
  \code{\link{biomarkerPowerTable}}
}
\note{See \code{\link{biomarkerPowerTable}} for examples.}
\keyword{design}
