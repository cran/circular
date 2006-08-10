\name{mle.wrappednormal}
\title{Wrapped Normal Maximum Likelihood Estimates}
\alias{mle.wrappednormal}
\alias{print.mle.wrappednormal}

\description{
Computes the maximum likelihood estimates for the parameters of a
Wrapped Normal distribution:  mean and concentration parameter.
}

\usage{
mle.wrappednormal(x, mu = NULL, rho = NULL, sd = NULL, K = NULL, tol = 1e-05, 
    min.sd = 1e-3, min.k = 10, max.iter = 100, verbose = FALSE, control.circular=list())
\method{print}{mle.wrappednormal}(x, digits = max(3, getOption("digits") - 3), \dots)
}

\arguments{
  \item{x}{a vector. The object is coerced to class
    \code{\link{circular}}.}
  \item{mu}{if \code{NULL} the maximum likelihood estimate of the mean
    direction is calculated, otherwise the value is coerced to an object of class \code{circular}.}
  \item{rho}{if {NULL} the maximum likelihood estimate of the
    concentration parameter is calculated.}
  \item{sd}{standard deviation of the (unwrapped) normal. Used as an
    alternative parametrization.}
  \item{K}{number of terms to be used in approximating the density.}
  \item{tol}{precision of the estimation.}
  \item{min.sd}{minimum value should be reached by the search procedure
    for the standard deviation parameter.}
  \item{min.k}{minimum number of terms used in approximating the density.}
  \item{max.iter}{maximum number of iterations.}
  \item{verbose}{logical, if \code{TRUE} information on the convergence
    process are printed.}
  \item{control.circular}{the attribute of the resulting objects (\code{mu})}
  \item{digits}{integer indicating the precision to be used.}
  \item{\dots}{further arguments passed to or from other methods.}  
}

\value{
  Returns a list with the following components:

  \item{call}{the match.call().}
  \item{mu}{the estimate of the mean direction or the value supplied as an object of class \code{circular}.}
  \item{rho}{the estimate of the concentration parameter or the
    value supplied}
  \item{sd}{the estimate of the standard deviation or the value supplied.}
  \item{est.mu}{TRUE if the estimator is reported.}
  \item{est.rho}{TRUE if the estimator is reported.}
  \item{convergence}{TRUE if the convergence is achieved.}
}

\author{Claudio Agostinelli}

\references{
  Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular
  Statistics, Section 4.2.1, World Scientific Press, Singapore.
}

\seealso{
\code{\link{mean.circular}}
}

\examples{
x <- rwrappednormal(n=50, mu=circular(0), rho=0.5)
mle.wrappednormal(x) # estimation of mu and rho (and sd)
mle.wrappednormal(x, mu=circular(0)) # estimation of rho (and sd) only
}

\keyword{htest}
