\name{var}
\title{Variance}
\alias{var}
\alias{var.default}
\alias{var.data.frame}

\description{
The \code{var} function from the \pkg{stats} is replace by a new
\code{method} in order to report the variance of circular data
appropriatly. \code{var.default} is an alias of the original function
\code{var} see \code{\link[stats]{cor}}. The behavior would be the same for objects which are not
from \code{\link{class}} \code{\link{data.frame}} and
\code{\link{circular}} (in the last case the variance is define as one
minus the mean resultant length divided by the sample size of data, see
\code{\link{var.circular}} for more details). The method for
\code{data.frame} will apply the \code{var} function to each columns.
}

\usage{
var(x, \dots)
\method{var}{default}(x, y = NULL, na.rm = FALSE, use, \dots)
\method{var}{data.frame}(x, \dots)
}

\arguments{
  \item{x}{a numeric vector, matrix or data frame.}
  \item{y}{\code{NULL} (default) or a vector, matrix or data frame with
    compatible dimensions to \code{x}.  The default is equivalent to
    \code{y = x} (but more efficient).}
  \item{na.rm}{logical. Should missing values be removed?}
  \item{use}{an optional character string giving a
    method for computing covariances in the presence
    of missing values.  This must be (an abbreviation of) one of the strings
    \code{"all.obs"}, \code{"complete.obs"} or \code{"pairwise.complete.obs"}.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\seealso{
\code{\link[stats]{cor}}, \code{\link{var.circular}}, \code{\link{rho.circular}} and \code{\link{summary.circular}}.
}

\keyword{univar}
