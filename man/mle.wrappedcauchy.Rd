\name{mle.wrappedcauchy}
\title{Wrapped Cauchy Maximum Likelihood Estimates}
\alias{mle.wrappedcauchy}
\alias{print.mle.wrappedcauchy}

\description{
Computes the maximum likelihood estimates for the parameters of a
Wrapped Cauchy distribution:  mean and concentration parameter.
}

\usage{
mle.wrappedcauchy(x, mu, rho, tol = 1e-15, max.iter = 100)
\method{print}{mle.wrappedcauchy}(x,
        digits = max(3, getOption("digits") - 3), \dots)
}

\arguments{
  \item{x}{a vector. The object is coerced to class
    \code{\link{circular}}.}
  \item{mu}{if missing the maximum likelihood estimate of the mean
    direction is calculated.}
  \item{rho}{if missing the maximum likelihood estimate of the
    concentration parameter is calculated.}
  \item{tol}{precision of the estimation.}
  \item{max.iter}{maximum number of iterations.}
   \item{digits}{integer indicating the precision to be used.}
  \item{\dots}{further arguments passed to or from other methods.}
 }

\value{
  Returns a list with the following components:

  \item{call}{the match.call().}
  \item{mu}{the estimate of the mean direction or the value supplied.}
  \item{rho}{the estimate of the concentration parameter or the
    value supplied}
  \item{convergence}{TRUE if convergence is achieved.}
}

\author{Claudio Agostinelli and Ulric Lund}

\references{
  Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular
  Statistics, Section 4.2.1, World Scientific Press, Singapore.
}

\seealso{
\code{\link{mean.circular}}
}

\examples{
x <- rwrappednormal(n=50, mu=0, rho=0.5)
mle.wrappednormal(x) # estimation of mu and rho
mle.wrappednormal(x, mu=0) # estimation of rho only
}

\keyword{htest}
