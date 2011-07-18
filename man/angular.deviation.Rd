\name{angular.deviation}
\title{A measure of deviation for Circular Data}
\alias{angular.deviation}

\description{
Returns the square root of twice one minus the mean resultant length divided by the sample size of a vector of circular data.
}

\usage{
angular.deviation(x, na.rm = FALSE)
}

\arguments{
  \item{x}{a vector. The object is coerced to class
    \code{\link{circular}}.}
  \item{na.rm}{logical, indicating if \code{\link{NA}}'s should be
    omitted.}
}

\value{
Returns the square root of twice one minus the mean resultant length divided by the sample size.
}

\author{Claudio Agostinelli}

\references{
Batchelet, E. (1981) Circular Statistics in Biology. Academic Press, London.
 
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 1.3, World Scientific Press, Singapore.

Zar, J.H. (2010) Biostatistical Analysis. Fifth edition. Pearson
Educational International.}

\seealso{
\code{\link{sd.circular}}, \code{\link{angular.variance}}, \code{\link{mean.circular}}, \code{\link{rho.circular}} and \code{\link{summary.circular}}.
}

\examples{
x <- rvonmises(n=100, mu=circular(0), kappa=1)
angular.deviation(x)
}

\keyword{univar}
