\name{median.circular}
\alias{median.circular}
\title{
Median Direction
}
\description{
Sample median direction for a vector of circular data
}
\usage{
\method{median}{circular}(x, na.rm = FALSE, ...)
}
\arguments{
  \item{x}{a vector. The object is coerced to class
    \code{\link{circular}}.
  }
  \item{na.rm}{logical, indicating if \code{\link{NA}}'s should be
	omitted.
  }
  \item{\dots}{NotYetUsed.}
}
\details{
The Definition in equations 2.32 & 2.33 from N.I. Fisher's 'Statistical Analysis of Circular Data', Cambridge Univ. Press 1993. is implemented.
Since version 0.4-4, the algorithm (not the definition) for the
calculation of the median is changed. For a measure of spread associated to the circular median use function \code{\link{meandeviation}}.
}

\value{
A scalar with the circular median value.

The median is returned as an object of class \code{circular}.
}

\references{
  N.I. Fisher (1993) Statistical Analysis of Circular Data, Cambridge
  University Press.
  
  R.Y. Liu and K. Singh (1992) Ordering Directional Data: Concepts of
  Data Depth on Circles and Spheres, The Annals of Statistics, vol. 20,
  n. 3, 1468-1484.
}

\author{
Claudio Agostinelli and Alessandro Gagliardi
}

\seealso{
\code{\link{meandeviation}}, \code{\link{mean.circular}}, \code{\link{var.circular}}, \code{\link{summary.circular}}, \code{\link{rho.circular}} and  \code{\link{medianHL.circular}}.
}
\examples{
# Compute the median direction of a random sample of observations.
x <- circular(runif(50, circular(0), pi))
median(x) #only the median is returned
meandeviation(x) #mean deviation is reported
}
\keyword{univar}
