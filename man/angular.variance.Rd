\name{angular.variance}
\title{A measure of variance for Circular Data}
\alias{angular.variance}

\description{
Returns twice one minus the mean resultant length divided by the sample size of a vector of circular data.
}

\usage{
angular.variance(x, na.rm = FALSE)
}

\arguments{
  \item{x}{a vector. The object is coerced to class
    \code{\link{circular}}.}
  \item{na.rm}{logical, indicating if \code{\link{NA}}'s should be
    omitted.}
}

\value{
Returns twice one minus the mean resultant length divided by the sample size.
}

\author{Claudio Agostinelli}

\references{
Batchelet, E. (1981) Circular Statistics in Biology. Academic Press, London.
 
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 1.3, World Scientific Press, Singapore.

Zar, J.H. (2010) Biostatistical Analysis. Fifth edition. Pearson
Educational International.
}

\seealso{
\code{\link{var.circular}}, \code{\link{angular.deviation}}, \code{\link{mean.circular}}, \code{\link{rho.circular}} and \code{\link{summary.circular}}.
}

\examples{
x <- rvonmises(n=100, mu=circular(0), kappa=1)
angular.variance(x)
}

\keyword{univar}
