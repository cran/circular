\name{range.circular}
\title{Circular Range}
\alias{range.circular}

\description{
Computes the circular range of a data set and performs a test of uniformity if specified.
}

\usage{
\method{range}{circular}(x, test=FALSE, na.rm = FALSE, finite = FALSE, 
  control.circular=list(), ...)
}

\arguments{
\item{x}{a vector. The object is coerced to class \code{\link{circular}}.}
\item{test}{logical flag:  if TRUE then the test of uniformity is
  performed; otherwise the test is not performed.  Default is FALSE.}
\item{na.rm}{logical, indicating if \code{\link{NA}}'s should be omitted.}
\item{finite}{logical, indicating if all non-finite elements should be
  omitted.}
\item{control.circular}{the attribute of the resulting object.}
\item{\ldots}{further parameter passed from/to the method.}
}

\value{
Returns the circular range as a \code{circular} object.  If the significance test is requested the p-value of the test is returned as p.value.
}

\details{
The circular range is the shortest arc on the circle containing the
entire set of data.  The p-value is computed using the exact
distribution of the circular range under the hypothesis of uniformity,
details can be found in Mardia and Jupp (1999) pag. 107.
}

\author{Claudio Agostinelli and Ulric Lund}

\references{
  K.V. Mardia and P.E. Jupp (1999) Directional Statistics, Wiley.
}

\seealso{
\code{\link{kuiper.test}}, \code{\link{rao.spacing.test}},
\code{\link{rayleigh.test}} and \code{\link{watson.test}}.
}

\examples{
data <- rvonmises(n=50, mu=circular(0), kappa=2)
range(data, test=TRUE)
data <- circular(runif(50, 0, 2*pi))
range(data, test=TRUE)
}

\keyword{univar}
\keyword{htest}
