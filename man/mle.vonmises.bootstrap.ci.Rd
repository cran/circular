\name{mle.vonmises.bootstrap.ci}
\title{Bootstrap Confidence Intervals}
\alias{mle.vonmises.bootstrap.ci}
\alias{print.mle.vonmises.bootstrap.ci}

\description{
Generates simple bootstrap confidence intervals for the parameters of a von Mises distribtution: the mean direction mu, and the concentration parameter kappa.
}


\usage{
mle.vonmises.bootstrap.ci(x, mu = NULL, bias = FALSE, alpha = 0.05, 
          reps = 1000, control.circular = list())
\method{print}{mle.vonmises.bootstrap.ci}(x, \dots)
}

\arguments{
  \item{x}{vector of angular measurements as a \code{circular} object.}
  \item{mu}{If \code{NULL} the value is estimated. This value is used in the bootstrap replications for \code{kappa}.}
\item{bias}{logical, if \code{TRUE}, the replication estimates for kappa are computed with a bias corrected method.  See \code{\link{mle.vonmises}}.  Default is \code{FALSE}, i.e. no bias correction.}
\item{alpha}{parameter determining level of confidence intervals. 1-alpha confidence intervals for \code{mu} and \code{kappa} are computed.  By default, 95\% confidence intervals are generated.}
\item{reps}{number of resampled data sets to use.  Default is 1000.}
  \item{control.circular}{the attribute of the resulting objects (\code{mu}, \code{mu.ci}).}
\item{\dots}{arguments passed to \code{\link{print.default}}.}
}

\value{
  A list is returned with the following components:

  \item{mu.ci}{limits of the confidence interval for mu as a \code{circular} object.}
  \item{kappa.ci}{limits of the confidence interval for kappa.}
  \item{mu}{estimate of mu as a \code{circular} object.}
  \item{kappa}{estimate of kappa.}
}

\details{
Percentile confidence intervals are computed by resampling from the original data set \code{reps} times.  For each resampled data set, the MLE's of mu and kappa are computed.  The bootstrap confidence intervals are the alpha/2 and 1-alpha/2 percentiles of the \code{reps} MLE's computed for each resampled data set.
}

\author{Claudio Agostinelli and Ulric Lund}

\seealso{
\code{\link{mle.vonmises}}
}

\examples{
x <- rvonmises(n=25, mu=circular(0), kappa=3)
x.bs <- mle.vonmises.bootstrap.ci(x, alpha=.10)
par(mfcol=c(1,2))
rose.diag(x.bs$mu, bins=30, main=expression(mu))
hist(x.bs$kappa, main=expression(kappa))
}

\keyword{htest}
