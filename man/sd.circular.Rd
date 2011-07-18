\name{sd.circular}
\title{Circular Standard Deviation}
\alias{sd.circular}

\description{
Returns the circular standard deviation of a vector of circular data which is
defined as the square root of minus 2 times the log of the mean resultant
length divided by the number of observations.
}

\usage{
\method{sd}{circular}(x, na.rm = FALSE, \dots)
}

\arguments{
  \item{x}{a vector. The object is coerced to class
    \code{\link{circular}}.}
  \item{na.rm}{logical, indicating if \code{\link{NA}}'s should be
    omitted.}
%%%  \item{control.circular}{the attributes of the resulting object.}  
  \item{\dots}{further arguments passed to or from other methods.}
}

\details{
Computes the circular standard deviation as defined by Mardia (1972) \deqn{\sqrt{-2\ln r}}{sqrt(-2 * ln(r))} where \code{r} is the mean resultant length of the data.
}

\value{
  Returns the circular standard deviation.

 %%%%% of the data as an object of class \code{\link{circular}} with the attributes given by \code{control.circular} or taken from \code{x} if missing in \code{control.circular}.
}

\author{Claudio Agostinelli and Jean-Olivier Irisson}

\references{
Mardia, K.V. (1972) Statistics of Directional Data. Academic Press,
London, sec. 26.5, p. 617

Fisher, N.I. (1993) Statistical analysis of circular data. Cambridge University Press.
  
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 1.3, World Scientific Press, Singapore.

Zar, J H (2010). Biostatistical analysis. Prentice Hall. sec. 26.5, p. 617
}

\seealso{
\code{\link{var.circular}}, \code{\link{angular.deviation}}, \code{\link{mean.circular}}, \code{\link{rho.circular}} and \code{\link{summary.circular}}.
}

\examples{
# Compute the circular standard deviation of a random
# sample of observations from a von Mises distribution
x <- rvonmises(n=100, mu=circular(0, units="degrees"), kappa=10)
sd(x)
}

\keyword{univar}
