\name{kuiper.test}
\title{Kuiper's Test}
\alias{kuiper.test}
\alias{print.kuiper.test}
\description{
Performs Kuiper's one-sample test of uniformity on the circle. 
}

\usage{
kuiper.test(x, alpha=0)
\method{print}{kuiper.test}(x, digits = 4, \dots)
}

\arguments{
  \item{x}{a vector. The object is coerced to class
    \code{\link{circular}}.}

  \item{alpha}{significance level of the test.  Possible levels are
    0.15, 0.1, 0.05, 0.025, 0.01.  Alpha may be omitted or set to zero,
    in which case a range for the p-value of the test will be printed.}
  \item{digits}{integer indicating the precision to be used.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{
A list with the statistic and alpha value.
}

\note{
Kuiper's one-sample test of uniformity is performed, and the results are printed to the screen.  If alpha is specified and non-zero, the test statistic is printed along with the critical value and decision.  If alpha is omitted, the test statistic is printed and a range for the p-value of the test is given.
}

\details{
Kuiper's test statistic is a rotation-invariant Kolmogorov-type test statistic.  The critical values of a modified Kuiper's test statistic are used according to the tabulation given in Stephens (1970). 
}

\author{Claudio Agostinelli and Ulric Lund}

\references{
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 7.2, World Scientific Press, Singapore.

Stephens, M. (1970).  Use of the Kolmogorov-Smirnov, Cramer-von Mises and related statistics without extensive tables.  Journal of the Royal Statistical Society, B32, 115-122.
}

\seealso{
\code{\link{range.circular}}, \code{\link{rao.spacing.test}},
\code{\link{rayleigh.test}} and \code{\link{watson.test}}
}

\examples{
# Generate data from the uniform distribution on the circle.
data <- circular(runif(100, 0, 2*pi))
kuiper.test(data)
# Generate data from the von Mises distribution.
data <- rvonmises(n=100, mu=circular(0), kappa=3)
kuiper.test(data, alpha=0.01)
}

\keyword{htest}
