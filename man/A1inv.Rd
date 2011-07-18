\name{A1inv}
\title{Inverse of A1}
\alias{A1inv}

\description{
Inverse function of the ratio of the first and zeroth order Bessel functions of the first kind.  This function is used to compute the maximum likelihood estimate of the concentration parameter of a von Mises distribution.
}

\usage{
A1inv(x)
}

\arguments{
\item{x}{numeric value in the interval between 0 and 1.}
}

\value{
Returns the value k, such that A1inv(x) = k, i.e. A1(k) = x.
}

\details{
A1inv(0) = 0 and A1inv(1) = Inf.  This function is useful in estimating
the concentration parameter of data from a von Mises distribution. Our
function use the results in Best and Fisher (1981). Tables use tabulated
values by Gumbel, Greenwood and Durand (1953).
}

\author{Claudio Agostinelli}

\references{
     BEST, D.J. and FISHER, N.I. 1981. The bias of the maximum 
       likelihood estimators for the von Mises-Fisher concentration 
       parameters. Communications in Statistics, 10, 493-502.

     GUMBEL, E.J., GREENWOOD, J.A. AND DURAND, D.  1953.  The circular 
       normal distribution: theory and tables.  J. Amer. Statis. Assoc., 
       48, 131-152.
}

\seealso{
\code{\link{mle.vonmises}}, \code{\link{A1}}, \code{\link[base]{besselI}}.
}

\examples{
#Generate data from a von Mises distribution
data <- rvonmises(n=50, mu=circular(pi), kappa=4)
#Estimate the concentration parameter
s <- sum(sin(data))
c <- sum(cos(data))
mean.dir <- atan2(s, c)
kappa <- A1inv(mean(cos(data - mean.dir)))
}

\keyword{math}
