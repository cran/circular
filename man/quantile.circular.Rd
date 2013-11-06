\name{quantile.circular}
\title{Sample Circular Quantiles}
\alias{quantile.circular}
\description{
  The function \code{quantile.circular} produces sample circular quantiles
  corresponding to the given probabilities for a circular data set.
}
\usage{
\method{quantile}{circular}(x, probs = seq(0, 1, 0.25), na.rm=FALSE, names = TRUE, type = 7, ...)
}
\arguments{
  \item{x}{numeric circular vector whose sample quantiles are wanted. \code{\link{NA}} and \code{NaN} values are not
    allowed in numeric vectors unless \code{na.rm} is \code{TRUE}.}
  \item{probs}{numeric vector of probabilities with values in
    \eqn{[0,1]}.  (Values up to \samp{2e-14} outside that
    range are accepted and moved to the nearby endpoint.)}
  \item{na.rm}{logical; if true, any \code{\link{NA}} and \code{NaN}'s
    are removed from \code{x} before the quantiles are computed.}
  \item{names}{logical; if true, the result has a \code{\link{names}}
    attribute.  Set to \code{FALSE} for speedup with many \code{probs}.}
  \item{type}{an integer between 1 and 9 selecting one of the
    nine quantile algorithms detailed below to be used.}
  \item{...}{further arguments passed to or from other methods.
		Like \code{quantile} and so on.}
}
\details{
  A vector of length \code{length(probs)} is returned;
  if \code{names = TRUE}, it has a \code{\link{names}} attribute.

  \code{\link{NA}} and \code{\link{NaN}} values in \code{probs} are
  propagated to the result.

  The algorithm will proceed how described below:
   1) Linearize the circular observations.
   2) Calulate the linear median like type establish.
   3) The value it will transformed in circular.
}

\section{Types}{
  See description on documentation of \code{quantile}.
}
\author{
  Claudio Agostinelli and Alessandro Gagliardi.
}

\examples{
x <- rvonmises(1001, mu=circular(pi), kappa=5)
quantile.circular(x) # Extremes & Quartiles by default
}
\keyword{univar}
