\name{totalvariation.circular}
\alias{totalvariation.circular}
\title{
   Conditional total variation distance between two circular samples. 
}
\description{
    The total variation distance between two circular samples is
    evaluated conditional on a circular modal region.    
}
\usage{
totalvariation.circular(x, y, z = NULL, q = 0.95, bw, adjust = 1,
type = c("K", "L"), kernel = c("vonmises", "wrappednormal"),
na.rm = FALSE, step = 0.001, eps.lower = 10^(-4), eps.upper = 10^(-4), ...)
}
\arguments{
  \item{x}{numeric or an object of class \code{\link{circular}}.}
  \item{y}{numeric or an object of class \code{\link{circular}}.}
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
  \item{\dots}{further arguments passed to the
    \code{modal.region.circular} function. Not used at present.}
}

\value{
A list of class \code{totalvariation.circular} with the following
components

\item{tv}{the (conditional) total variation.}
\item{ovl}{the (conditional) overlapping coefficient.}
\item{q}{the order of the modal regions.}
\item{bw}{the bandwidth value as in input.} 
\item{modal.x}{an obejct of class \code{\link{modal.region.circular}}
  for the \code{x} data set.}
\item{modal.y}{an obejct of class \code{\link{modal.region.circular}}
  for the \code{y} data set.}
\item{density.x}{an obejct of class \code{\link{density.circular}}
  for the \code{x} data set.}
\item{density.y}{an obejct of class \code{\link{density.circular}}
  for the \code{y} data set.}
\item{density}{a function which report the positive part of the
  difference between the estimated density of the two data sets.}
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
\code{\link{modal.region.circular}}
}
\examples{
x <- rvonmises(100, circular(pi), 10)
y <- rvonmises(100, circular(pi+pi/8), 10)
res <- totalvariation.circular(x,y,bw=50)
plot(res)
}

\keyword{univariate}
