\name{medianCircular}
\alias{medianCircular}
\title{
Median Direction
}
\description{
Sample median direction and corresponding deviation measure for a vector of circular data
}
\usage{
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
For now only the definition in equations 2.32 & 2.33 from N.I. Fisher's 'Statistical Analysis of Circular Data', Cambridge Univ. Press 1993. is implemented.
}
\value{
If \code{dispersion==FALSE} a scalar with the median of the data otherwise a list with the following two components
\item{median}{the median.}
\item{deviation}{the associated deviation from the median.} 

The median is returned as an object of class \code{circular} with the attribute given by \code{control.circular} or from \code{x} if missed in \code{control.circular}.
}

\references{
  N.I. Fisher (1993) Statistical Analysis of Circular Data, Cambridge University Press.
}
\author{
Claudio Agostinelli
}

\seealso{
\code{\link{mean.circular}}, \code{\link{var.circular}}, \code{\link{summary.circular}} and \code{\link{rho.circular}}.
}
\examples{
# Compute the median direction of a random sample of observations.
x <- circular(runif(50, circular(0), pi))
medianCircular(x) #only the median is returned
medianCircular(x, deviation=TRUE) #both median and deviation are reported
}
\keyword{univar}
