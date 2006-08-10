\name{rayleigh.test}
\title{Rayleigh Test of Uniformity}
\alias{rayleigh.test}
\alias{print.rayleigh.test}

\description{
Performs a Rayleigh test of uniformity, assessing the significance of
the mean resultant length.  The alternative hypothesis is a unimodal
distribution with unknown mean direction and unknown mean resultant
length if \code{mu} is \code{NULL} otherwise the alternative hypothesis is a unimodal distribution with a specified mean direction and unknown mean resultant length.  
}

\usage{
rayleigh.test(x, mu = NULL)
\method{print}{rayleigh.test}(x, digits=4, \dots)
}

\arguments{
  \item{x}{a vector. The object is coerced to class
    \code{\link{circular}}.}
  \item{mu}{Specified mean direction in alternative hypothesis as a \code{circular} object.}
    \item{digits}{integer indicating the precision to be used.}
	\item{\dots}{further arguments passed to or from other methods.}
  }

\value{
Returns a list with three components:  the mean resultant length, \code{statistic}, the p-value of the test statistic, \code{p.value} and
the value of the alternative mean direction \code{mu}.
}

\author{Claudio Agostinelli and Ulric Lund}

\references{
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Sections 3.3.2 and 3.4.1, World Scientific Press, Singapore.
}

\seealso{
\code{\link{range.circular}}, \code{\link{kuiper.test}},
\code{\link{rao.spacing.test}} and \code{\link{watson.test}}
}


\examples{
x <- rvonmises(n=25, mu=circular(pi), kappa=2)
# General alternative
rayleigh.test(x)
# Specified alternative
rayleigh.test(x, mu=circular(0))
}

\keyword{htest}
