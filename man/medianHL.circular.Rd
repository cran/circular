\name{medianHL.circular}
\alias{medianHL.circular}
\alias{medianHL}
\alias{medianHL.default}
\title{
  Median using Hodges-Lehmann estimate.
}
\description{
  Sample median for a vector of data using Hodges-Lehmann estimate
  and Sample median direction measure for a vector of circular data
  using Hodges-Lehmann estimate.
}
\usage{
medianHL(x, na.rm=FALSE, ...)
\method{medianHL}{default}(x, na.rm=FALSE,
       method=c("HL1","HL2","HL3"), prop=NULL)
\method{medianHL}{circular}(x, na.rm=FALSE,
       method=c("HL1","HL2","HL3"), prop=NULL)
}
\arguments{
  \item{x}{a vector. For the function \code{medianHL.circular} the object is
    coerced to class \code{\link{circular}}.}
  \item{na.rm}{logical, indicating if \code{\link{NA}}'s should be omitted.}
  \item{method}{The method used to calculate the median, see details below.}
  \item{prop}{The proportion of pairs that are sampled. If \code{NULL}
    all combinations are used. It must be a number in the interval (0,1)
    or \code{NULL}.}
  \item{...}{further arguments passed to the next method.}
}
\details{
The algorithm is as follows:

The algorithm will create pairs of elements of the vector \code{x}.

It will calculate the circular mean on those pairs.

It will calculate the circular median on these averages.

The type of pairs considered are controlled by \code{method}:

if \code{method} is "HL1" are considered unordered pairs without
replications and repetition in the number of \code{(n*(n-1))/2} pairs;

if \code{method} is "HL2" are considered unordered pairs without
replications in the number of \code{(n*(n+1))/2} pairs;

if \code{method} is "HL3" all pairs are considered in the number of \code{n^2}.

If \code{prop} is not \code{NULL}, the algorithm will consider a
subsample following the rules specified by \code{method}, however, the
number of pairs considered is prop * (number of pairs defined by \code{method}).

For more details see Bennett Sango Otieno, 'An Alternative Estimate of
Preferred Direction for Circular Data', Virginia Tech (2002) pag. 27-28 and 46-47.
}

\value{
  For \code{medianHL.circular} the median is returned as an object of
  class \code{circular} with the attribute given by those of \code{x}.
  An attributes \code{medians} reports all the averages which are
  minimizer of the circular median function.   
}

\references{
  Bennett Sango Otieno, An Alternative Estimate of Preferred Direction for Circular Data, Virginia Tech (July 2002).

  Bennett Sango Otieno and Christine M. Anderson-Cook,Measures of preferred direction for environmental and ecological circular data, Springer (June 2004).
}

\author{
  Claudio Agostinelli and Alessandro Gagliardi.
}

\seealso{
  \code{\link{mean.circular}}, \code{\link{median.circular}}.
}
\examples{
# Compute the median direction of a random sample of observations.
  x <- circular(runif(50, circular(0), pi))
# Calculate the three medians for each method without \code{prop} argument.
  medianHL.circular(x,method="HL1")
  medianHL.circular(x,method="HL2")
  medianHL.circular(x,method="HL3")
}
\keyword{univar}
