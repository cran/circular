\name{rao.spacing.test}
\title{Rao's Spacing Test of Uniformity}
\alias{rao.spacing.test}
\alias{print.rao.spacing.test}

\description{
Performs Rao's spacing test of uniformity. 
}

\usage{
rao.spacing.test(x, alpha=0)
\method{print}{rao.spacing.test}(x, digits = 4, \dots)
}

\arguments{
  \item{x}{a vector. The object is coerced to class
    \code{\link{circular}}.}
  \item{alpha}{numeric value specifying the significance level of the
  test. The default value is 0, in which case, a range for the p-value
  will be returned.  Valid significance levels are 0.10, 0.05, 0.01 and
  0.001.}
\item{digits}{integer indicating the precision to be used.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{
a list with the statistic, alpha and the number of observations.
}

\details{
If alpha is specified, critical values are determined from a table of simulated critical points (see reference below).  If alpha is not specified, a range for the p-value is determined using the table of simulated critical points.
}

\author{Claudio Agostinelli and Ulric Lund}

\references{
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 7.4, World Scientific Press, Singapore.

Rao, J.S. (1976). Some tests based on arc-lengths for the circle. Sankhya, The Indian Journal of Statistics, Serial B(4), 38, 329–338.

Russell, G.S. and Levitin, D.J. (1995).  An expanded table of probability values for Rao's Spacing Test.  Communications in Statistics - Simulation and Computation, 24, 4, 879–888.
}

\seealso{
\code{\link{range.circular}}, \code{\link{kuiper.test}}, \code{\link{rayleigh.test}} and \code{\link{watson.test}}
}

\keyword{htest}
