\name{var.circular}
\title{A measure of variance for Circular Data}
\alias{var.circular}

\description{
Returns one minus the mean resultant length divided by the sample size of a vector of circular data. 
}

\usage{
\method{var}{circular}(x, na.rm = FALSE, \dots)
}

\arguments{
  \item{x}{a vector. The object is coerced to class
    \code{\link{circular}}.}
  \item{na.rm}{logical, indicating if \code{\link{NA}}'s should be
    omitted.}
  \item{\dots}{further arguments passed to or from other methods.}
}


\value{
Returns one minus the mean resultant length divided by the sample size of data.
}

\author{Claudio Agostinelli and Ulric Lund}

\references{
Mardia, K.V. (1972) Statistics of Directional Data. Academic Press, London. 

Fisher, N.I. (1993) Statistical analysis of circular data. Cambridge
University Press.

Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 1.3, World Scientific Press, Singapore.
}

\seealso{
\code{\link{sd.circular}}, \code{\link{angular.variance}}, \code{\link{mean.circular}}, \code{\link{rho.circular}} and \code{\link{summary.circular}}.
}

\examples{
x <- rvonmises(n=100, mu=circular(0), kappa=1)
var(x)
}

\keyword{univar}
