\name{meandeviation}
\title{A measure of deviation for Circular Data}
\alias{meandeviation}

\description{
Returns a measure of spread associated with the circular median of a vector of circular data.
}

\usage{
meandeviation(x, na.rm = FALSE)
}

\arguments{
  \item{x}{a vector. The object is coerced to class
    \code{\link{circular}}.}
  \item{na.rm}{logical, indicating if \code{\link{NA}}'s should be
    omitted.}
}

\details{
See equation (2.33) at pag. 36 in Fisher (1993)
for its definition. In the case the circular median is not defined, that
is, every angle is a median axis, the mean deviation is not reported and
set to \code{NA}.
}

\value{
Returns a measure of spread associated with the circular median of a
vector of circular data.
}

\author{Claudio Agostinelli and Alessandro Gagliardi}

\references{
  N.I. Fisher (1993) Statistical Analysis of Circular Data, Cambridge
  University Press.
}

\seealso{
\code{\link{median.circular}}, \code{\link{sd.circular}}, \code{\link{angular.variance}}, \code{\link{angular.deviation}}, \code{\link{mean.circular}}, \code{\link{rho.circular}} and \code{\link{summary.circular}}.
}

\examples{
x <- rvonmises(n=100, mu=circular(0), kappa=1)
meandeviation(x)
}

\keyword{univar}
