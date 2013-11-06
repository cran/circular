\name{medianCircular}
\alias{median.circular}
\alias{medianCircular}
\title{
Median Direction
}
\description{
Sample median direction and corresponding deviation measure for a vector of circular data
}
\usage{
\method{median}{circular}(x, na.rm = FALSE)
medianCircular(x, na.rm = FALSE, type = "Fisher", 
  deviation = FALSE, control.circular = list(), \dots)
}
\arguments{
  \item{x}{a vector. The object is coerced to class
    \code{\link{circular}}.
  }
  \item{na.rm}{logical, indicating if \code{\link{NA}}'s should be
	omitted.
  }
  \item{type}{the only implemented method is \code{Fisher}.}
  \item{deviation}{
  if \code{TRUE} a list with two objects is returned, the first component be the median, the second component be the deviation associated with the median.
  }
  \item{control.circular}{the attribute of the resulting object.}
  \item{\dots}{NotYetUsed.}
}
\details{
The Definition in equations 2.32 & 2.33 from N.I. Fisher's 'Statistical Analysis of Circular Data', Cambridge Univ. Press 1993. is implemented.
Since version 0.4-4, the algorithm (not the definition) for the
calculation of the median is changed. The function \code{medianCircular}
is deprecated and it will be remove in the next release of the
package. For a measure of spread associated to the circular median use
function \code{\link{meandeviation}}.
}

\value{
If \code{dispersion==FALSE} a scalar with the median of the data otherwise a list with the following two components
\item{median}{the median.}
\item{deviation}{the associated deviation from the median.} 

The median is returned as an object of class \code{circular} with the attribute given by \code{control.circular} or from \code{x} if missed in \code{control.circular}.
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
\dontrun{
  medianCircular(x) #only the median is returned
  medianCircular(x, deviation=TRUE) #both median and deviation are reported
}
}
\keyword{univar}
