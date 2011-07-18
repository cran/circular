\name{mle.vonmises}
\title{von Mises Maximum Likelihood Estimates}
\alias{mle.vonmises}
\alias{print.mle.vonmises}

\description{
Computes the maximum likelihood estimates for the parameters of a von Mises distribution:  the mean direction and the concentration parameter.
}

\usage{
mle.vonmises(x, mu=NULL, kappa=NULL, bias=FALSE, control.circular=list())
\method{print}{mle.vonmises}(x,
    digits = max(3, getOption("digits") - 3), \dots)
}

\arguments{
  \item{x}{a vector. The object is coerced to class
    \code{\link{circular}}.}
  \item{mu}{if \code{NULL} the maximum likelihood estimate of the mean
    direction is calculated. If provided it is coerced to a class \code{circular}.}
  \item{kappa}{if \code{NULL} the maximum likelihood estimate of the
    concentration parameter is calculated.}
  \item{bias}{logical, if \code{TRUE}, the estimate for kappa is
    computed with a bias corrected method.  Default is \code{FALSE},
    i.e. no bias correction.}
  \item{control.circular}{the attribute of the resulting objects (\code{mu})}
  \item{digits}{integer indicating the precision to be used.}
  \item{\dots}{further arguments passed to or from other methods.}  
}

\value{
  Returns a list with the following components:

  \item{call}{the \code{\link[base]{match.call}} result.}
  \item{mu}{the estimate of the mean direction or the value supplied as an object of class \code{circular}.}
  \item{kappa}{the estimate of the concentration parameter or the
    value supplied}
  \item{se.mu}{the standard error for the estimate of the mean
    direction (0 if the value is supplied) in the same units of \code{mu}.}
  \item{se.kappa}{the standard error for the estimate of the
    concentration parameter (0 if the value is supplied).}
  \item{est.mu}{TRUE if the estimator is reported.}
  \item{est.kappa}{TRUE if the estimator is reported.}
}

\details{
Best and Fisher (1981) show that the MLE of kappa is seriously biased when both sample size and mean resultant length are small. They suggest a bias-corrected estimate for kappa when n < 16.
}

\author{Claudio Agostinelli and Ulric Lund}

\references{
  Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular
  Statistics, Section 4.2.1, World Scientific Press, Singapore.

  
  Best, D. and Fisher N. (1981).  The bias of the maximum likelihood
  estimators of the von Mises-Fisher concentration parameters.
  Communications in Statistics - Simulation and Computation,
  B10(5), 493-502.
}

\seealso{
\code{\link{mean.circular}} and \code{\link{mle.vonmises.bootstrap.ci}}
}

\examples{
x <- rvonmises(n=50, mu=circular(0), kappa=5)
mle.vonmises(x) # estimation of mu and kappa
mle.vonmises(x, mu=circular(0)) # estimation of kappa only
}

\keyword{htest}
