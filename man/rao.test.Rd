\name{rao.test}
\title{Rao's Tests for Homogeneity}
\alias{rao.test}
\alias{print.rao.test}

\description{
Performs Rao's test for homogeneity on k populations of angular data. 
}

\usage{
rao.test(\dots, alpha=0)
\method{print}{rao.test}(x, digits = 4, \dots)
}

\arguments{
  \item{\dots}{a sequence of \code{\link{circular}} for the
	\code{rao.test} and further arguments passed to or from other
	methods for the \code{print.rao.test} function.}
  \item{alpha}{numeric value specifying the significance level of the
  test. Default is 0, in which case p-values for the test statistic is
  printed.}
  \item{x}{an object from the \code{rao.test}.}
  \item{digits}{integer indicating the precision to be used.}
}

\value{
A list with the statistic and p.value for the mean and the dispersion
and the value of alpha. 
}

\note{
The test is performed, and the results are written to the screen. Test results are given for both the test of equality of polar vectors, and of dispersions.  If alpha is specified, the test statistic is printed, along with the level critical value.  If alpha is not specified, a p-value for the test is printed.  
}

\details{
Critical values and p-values are determined according to the chi-squared approximation of the test statistic.
}

\author{Claudio Agostinelli and Ulric Lund}

\references{
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 7.6.1, World Scientific Press, Singapore.

Rao, J.S. (1967). Large sample tests for the homogeneity of angular data, Sankhya, Ser, B., 28.
}

\examples{
x <- rvonmises(100, circular(0), kappa=10)
y <- rvonmises(100, circular(0), kappa=10)

rao.test(x, y)

}


\keyword{htest}
