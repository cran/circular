\name{var.circular}
\title{A measure of variance for Circular Data}
\alias{var.circular}

\description{
Returns one minus the mean resultant length divided by the sample size of a vector of circular data. 
}

\usage{
\method{var}{circular}(x, na.rm = FALSE, only.var = TRUE, \dots)
}

\arguments{
  \item{x}{a vector. The object is coerced to class
    \code{\link{circular}}.}
  \item{na.rm}{logical, indicating if \code{\link{NA}}'s should be
    omitted.}
  \item{only.var}{logical, if \code{FALSE} the mean result length is
    reported as well.}
  \item{\dots}{further arguments passed to or from other methods.}
}


\value{
Returns one minus the mean resultant length divided by the sample size of data if \code{only.var} is
\code{TRUE}, otherwise a vector with the following components: the number
of observations, the mean resultant length, the mean resultant length
divided by the sample size, and one minus the mean resultant length divided by the sample size.
}

\author{Claudio Agostinelli and Ulric Lund}

\references{
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 1.3, World Scientific Press, Singapore.
}

\seealso{
\code{\link{mean.circular}}, \code{\link{rho.circular}} and \code{\link{summary.circular}}.
}

\examples{
# Compute the mean resultant length of a random sample of observations.
x <- rvonmises(n=100, mu=circular(0), kappa=1)
var(x)
}

\keyword{univar}
