\name{rho.circular}
\title{Mean Resultant Length}
\alias{rho.circular}

\description{
Returns the mean resultant length of a vector of circular data. 
}

\usage{
rho.circular(x, na.rm = FALSE)
}

\arguments{
  \item{x}{a vector. The object is coerced to class
    \code{\link{circular}}.}
  \item{na.rm}{logical, indicating if \code{\link{NA}}'s should be omitted.}
}


\value{
Returns the mean resultant length of data.
}

\details{
Each observation is treated as a unit vector, or point  on the unit circle.  The resultant vector of the observations is found, and the length of the resultant vector divided by the sample size is returned.
}

\author{Claudio Agostinelli and Ulric Lund}

\references{
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 1.3, World Scientific Press, Singapore.
}

\seealso{
\code{\link{mean.circular}}, \code{\link{var.circular}},
\code{\link{summary.circular}} and \code{\link{mle.vonmises}}.
}

\examples{
# Compute the mean resultant length of a random sample of observations.
data <- circular(runif(100, 0, 2*pi))
rho.circular(data)
}

\keyword{univar}
