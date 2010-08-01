\name{dist.circular}
\alias{dist.circular}
\concept{dissimilarity}
\title{Distance Matrix Computation for Circular Data}
\description{
  This function computes and returns the distance matrix computed by
  using the specified distance measure to compute the distances between
  the rows of a data matrix containing circular data.
}
\usage{
dist.circular(x, method = "correlation", diag = FALSE, upper = FALSE)
}

\arguments{
  \item{x}{a numeric matrix of class \code{\link{circular}}.}
  \item{method}{the distance measure to be used. This must be one of
    \code{"correlation"}, \code{"angularseparation"}, \code{"chord"},
    \code{"geodesic"}. Any unambiguous substring can be given.}
  \item{diag}{logical value indicating whether the diagonal of the
    distance matrix should be printed by \code{print.dist}.}
  \item{upper}{logical value indicating whether the upper triangle of the
    distance matrix should be printed by \code{print.dist}.}
}

\details{
  Available distance measures are (written for two vectors \eqn{x} and
  \eqn{y}):
  \describe{
    \item{\code{correlation}:}{\eqn{\sqrt{1 - \rho}}{sqrt(1-rho)} where \eqn{\rho}{rho} is the Circular Correlation coefficient defined as \deqn{\frac{\sum_{i=1}^n \sin(x_i - \mu_x) \sin(y_i - \mu_y)}{\sqrt{\sum_{i=1}^n \sin^2(x_i - \mu_x) \sum_{i=1}^n \sin^2(y_i - \mu_y)}}}{sum(sin(x - mux)*sin(y - muy))/(sum(sin(x - mux)^2)*sum(sin(y - muy)^2))^(1/2)} and \eqn{\mu_x}{mux}, \eqn{\mu_y}{muy} are the mean direction of the two vectors} 
    \item{\code{angularseparation}:}{\eqn{\sum_{i=1}^n 1 - cos(x_i - y_i)}{sum(1 - cos(x - y))}}

    \item{\code{chord}:}{\eqn{\sum_{i=1}^n \sqrt{2 (1 - \cos(x_i - y_i))}}{sum(sqrt(2 (1 - cos(x - y))))}}
    \item{\code{geodesic}:}{\eqn{\sum_{i=1}^n \pi - |\pi - |x_i - y_i||}{sum(pi - abs(pi - abs(x - y)))} where the abs(x - y) is expressed with an angle in [-pi,pi]}
  }

  Missing values are allowed, and are excluded from all computations
  involving the rows within which they occur.
  Further, when \code{Inf} values are involved, all pairs of values are
  excluded when their contribution to the distance gave \code{NaN} or
  \code{NA}. \cr
  If some columns are excluded in calculating the sum is scaled up proportionally
  to the number of columns used.  If all pairs are excluded when calculating a
  particular distance, the value is \code{NA}.

}
\value{
  \code{dist.circular} returns an object of class \code{"dist"}.

  The lower triangle of the distance matrix stored by columns in a
  vector, say \code{do}. If \code{n} is the number of
  observations, i.e., \code{n <- attr(do, "Size")}, then
  for \eqn{i < j <= n}, the dissimilarity between (row) i and j is
  \code{do[n*(i-1) - i*(i-1)/2 + j-i]}.
  The length of the vector is \eqn{n*(n-1)/2}, i.e., of order \eqn{n^2}.

  The object has the following attributes (besides \code{"class"} equal
  to \code{"dist"}):
  \item{Size}{integer, the number of observations in the dataset.}
  \item{Labels}{optionally, contains the labels, if any, of the
    observations of the dataset.}
  \item{Diag, Upper}{logicals corresponding to the arguments \code{diag}
    and \code{upper} above, specifying how the object should be printed.}
  \item{call}{optionally, the \code{\link{call}} used to create the
    object.}
  \item{method}{optionally, the distance method used; resulting from
    \code{\link{dist.circular}()}, the (\code{\link{match.arg}()}ed) \code{method}
    argument.}
}

%\references{
%}

\seealso{
  \code{\link[stats]{dist}}
}
%\examples{
%
%}
\keyword{multivariate}
\keyword{cluster}
