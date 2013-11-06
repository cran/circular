\name{modal.region}
\alias{modal.region}
\alias{modal.region.default}
\alias{modal.region.circular}

\title{
   Modal regions
}
\description{
  Evaluate the modal regions for a data set. Only the version for
   circular data is implemented.
}
\usage{
modal.region(x, ...)
\method{modal.region}{default}(x, ...)
\method{modal.region}{circular}(x, z=NULL, q=0.95, bw, adjust = 1,
  type = c("K", "L"), kernel = c("vonmises", "wrappednormal"),
  na.rm = FALSE, step=0.01, eps.lower=10^(-4), eps.upper=10^(-4), ...)
}
\arguments{
  \item{x}{numeric or an object of class \code{\link{circular}}.}
  \item{z}{numeric or object of class \code{\link{circular}}. The grid
    were the kernel density estimate will be evaluated. If \code{NULL}
    equally spaced points in the interval [0,2*pi) with step \code{step}.}
  \item{q}{numeric in the interval [0,1]. The quantile of the modal
    region.} 
  \item{bw}{the smoothing bandwidth to be used. When the \code{kernel}
    is \code{vonmises} the bandwidth is equal to the concentration
    parameter.}
  \item{adjust}{the bandwidth used is actually \code{adjust*bw}. This
    makes it easy to specify values like ``half the default bandwidth''.}
  \item{type}{Not Yet Used.}
  \item{kernel}{a character string giving the smoothing kernel to be
          used. This must be one of \code{"vonmises"} or
      \code{"wrappednormal"}, that are kernels of \code{type} \code{"K"}.}
  \item{na.rm}{logical; if \code{TRUE}, missing values are removed from
    \code{x}. If \code{FALSE} any missing values cause an error.}
  \item{step}{numeric. Used in the construction of the regular grid \code{z}.}
  \item{eps.lower,eps.upper}{the cut point in the density is searched in
   the interval [min(density)*(1+eps.lower),max(density)*(1-eps.upper)].}
  \item{\dots}{further arguments passed to the next methods.}
}
\details{
  Only the version for circular data is actually implemented.
}
\value{
A list of class \code{modal.region.circular} with the following elements

\item{zeros}{extremes of modal regions, possible as a matrix}
\item{areas}{a list with two components: \code{tot} with the total
  (area under the density) probability, which should approximately equal
  to \code{q} and \code{areas} with the probability of each modal region.}
\item{density}{the object from function \code{density.circular}.}
\item{q}{the modal region order as in input.}
\item{level}{the cut point at the density scale.}
}

\references{
L.G.R. Oliveira-Santos, C.A. Zucco and C. Agostinelli (2013) Using
conditional circular kernel density functions to test hypotheses on
animal circadian activity. Animal Behaviour, 85(1) 269-280.
}

\author{
  Claudio Agostinelli
}

\seealso{
  \code{\link{totalvariation.circular}}
}
\examples{
x <- rvonmises(100, circular(pi), 10)
res <- modal.region(x, bw=50)
plot(res)
}
\keyword{univariate}
