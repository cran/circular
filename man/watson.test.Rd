\name{watson.test}
\title{Watson's Test }
\alias{watson.test}
\alias{print.watson.test}

\description{
Performs a Watson's goodness of fit test for the von Mises or circular uniform distribution. 
}

\usage{
watson.test(x, alpha=0, dist=c("uniform", "vonmises"))
\method{print}{watson.test}(x, digits = 4, \dots)
}

\arguments{
\item{x}{a vector. The object is coerced to class
    \code{\link{circular}}.}
\item{alpha}{significance level of the test.  Valid levels are 0.01, 0.05, 0.1.  This argument may be ommited, in which case, a range for the p-value will be returned.}
\item{dist}{distribution to test for.  The default is the uniform
  distribution. To test for the von Mises distribution, set \code{dist}
  to "vonmises".}
  \item{digits}{integer indicating the precision to be used.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{
a list with the statistic, alpha, the number of observations, the
distribution and 'row' which is used by \code{print.watson.test} to
evaluate the p-value.
}

\details{If \code{dist} = "uniform", Watson's one-sample test for the circular uniform distribution is performed, and the results are printed.  If alpha is specified and non-zero, the test statistic is printed along with the critical value and decision.  If alpha is omitted, the test statistic is printed and a range for the p-value of the test is given.

If \code{dist} = "vonmises", estimates of the population parameters are used to evaluate the von Mises distribution function at all data points, thereby arriving at a sample of approximately uniformly distributed data, if the original observations have a von Mises distribution. The one-sample Watson test is then applied to the transformed data as above.
}

\author{Claudio Agostinelli and Ulric Lund}

\references{
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 7.2, World Scientific Press, Singapore.

Stephens, M. (1970).  Use of the Kolmogorov-Smirnov, Cramer-von Mises and related statistics without extensive tables.  Journal of the Royal Statistical Society, B32, 115-122.
}

\seealso{
\code{\link{range.circular}}, \code{\link{kuiper.test}}, \code{\link{rao.spacing.test}} and \code{\link{rayleigh.test}}
}

\examples{
# Generate data from the uniform distribution on the circle.
x <- circular(runif(100, 0, 2*pi))  
watson.test(x)
# Generate data from a von Mises distribution.
x <- rvonmises(n=50, mu=circular(0), kappa=4) 
watson.test(x, alpha=0.05, dist="vonmises")
}

\keyword{htest}
