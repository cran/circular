\name{watson.two.test}
\title{Watson's Two-Sample Test of Homogeneity}
\alias{watson.two.test}
\alias{print.watson.two.test}

\description{
Performs Watson's test for homogeneity on two samples of circular data. 
}

\usage{
watson.two.test(x, y, alpha=0)
\method{print}{watson.two.test}(x, digits=4, \dots)
}

\arguments{
\item{x}{a vector. The object is coerced to class
    \code{\link{circular}}.}
\item{y}{a vector. The object is coerced to class
    \code{\link{circular}}.}
\item{alpha}{significance level of the test.  Valid levels are 0.001, 0.01, 0.05, 0.1.  This argument may be ommited, in which case, a range for the p-value will be returned.}
  \item{digits}{integer indicating the precision to be used.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{
a list with statistic, alpha and the number of observations of the first and 
second sample.
}

\details{
Watson's two-sample test of homogeneity is performed, and the results are printed.  If alpha is specified and non-zero, the test statistic is printed along with the critical value and decision.  If alpha is omitted, the test statistic is printed and a range for the p-value of the test is given.

Critical values for the test statistic are obtained using the asymptotic distribution of the test statistic.  It is recommended to use the obtained critical values and ranges for p-values only for combined sample sizes in excess of 17.  Tables are available for smaller sample sizes and can be found in Mardia (1972) for instance.
}

\author{Claudio Agostinelli and Ulric Lund}

\references{
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 7.5, World Scientific Press, Singapore.
}

\examples{
# Perform a two-sample test of homogeneity on two
# simulated data sets.
data1 <- rvonmises(n=20, mu=circular(0), kappa=3)
data2 <- rvonmises(n=20, mu=circular(pi), kappa=2)
watson.two.test(data1, data2, alpha=0.05)
watson.two.test(data1, data2)
}

\keyword{htest}
