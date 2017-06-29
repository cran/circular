\name{intersect.modal.region}
\alias{intersect.modal.region}
\alias{intersect.modal.region.default}
\alias{intersect.modal.region.circular}
  
\title{
  Intersection between model region and a given interval.
}
\description{
Find an estimates of the probability of the intersection between a modal region and a given interval. 
}
\usage{
intersect.modal.region(x, ...)
\method{intersect.modal.region}{default}(x, ...)
\method{intersect.modal.region}{circular}(x, breaks, z = NULL, q = 0.95, bw,
  adjust = 1, type = c("K", "L"), kernel = c("vonmises", "wrappednormal"),
  na.rm = FALSE, step = 0.01, eps.lower = 10^(-4), eps.upper = 10^(-4), ...)  
}
  
\arguments{
  \item{x}{numeric or an object of class \code{\link{circular}}.}
  \item{breaks}{a matrix with two columns. Each row specifies a sub-interval.}
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
For the circular method a list with the following three components  
\item{tot}{the total area.}
\item{areas}{information for each subinterval.}
\item{breaks}{the extremes of each subinterval.}
}
\author{
  Claudio Agostinelli
}

\seealso{
\code{\link{modal.region}}
}
\examples{
  x <- rvonmises(100, circular(pi), 10)  
  res <- intersect.modal.region(x, breaks=circular(matrix(c(pi,pi+pi/12,
    pi-pi/12, pi), ncol=2, byrow=TRUE)), bw=50)
  res$tot

  x <- rvonmises(100, circular(0), 10)
  res <- intersect.modal.region(x, breaks=circular(matrix(c(pi,pi+pi/12),
    ncol=2)), bw=50)
  res$tot
  
  res <- intersect.modal.region(x, breaks=circular(matrix(c(pi/12,
    2*pi-pi/12), ncol=2, byrow=TRUE)), bw=50)
  res$tot
}
