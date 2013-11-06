\name{summary.circular}
\title{Circular Summary Statistics}
\alias{summary.circular}

\description{
Computes circular summary statistics including the sample size, mean
direction and mean resultant length and quartiles.
}

\usage{
\method{summary}{circular}(object,
   digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{object}{an object of class \code{\link{circular}}.}
  \item{digits}{digits to be used in printing.}
  \item{...}{parameters passed to \code{\link{summary.matrix}} if needed.}
}

\value{
Returns a vector with the sample size, the sample mean direction and the sample mean resultant length.
}

\details{
Each observation is  treated  as  a  unit  vector or a point on the unit circle.  The resultant vector of the observations is found, and the direction of the resultant vector is returned as well as its length divided by the sample size.  
}

\author{Claudio Agostinelli, David Andel and Alessandro Gagliardi}

\references{
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 1.3, World Scientific Press, Singapore.
}

\seealso{
\code{\link{mean.circular}}, \code{\link{median.circular}},
\code{\link{quantile.circular}}, \code{\link{var.circular}},
\code{\link{mle.vonmises}}, \code{\link{rho.circular}}.
}

\examples{
# Compute summary statistics of a random sample of observations. 
data <- circular(runif(50, 0, pi))
summary(data)
summary(data.frame(data, runif(50, 0, pi)))
}

\keyword{univar}
