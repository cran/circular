\name{density.circular}
\alias{density.circular}
\alias{print.density.circular}

\title{Kernel Density Estimation for Circular Data}
\description{
  The function \code{density.circular} computes kernel density estimates
  with the given kernel and bandwidth for circular data.
}

\usage{
\method{density}{circular}(x, z=NULL, bw, adjust = 1, type = c("K", "L"),
     kernel = c("vonmises", "wrappednormal"), na.rm = FALSE, from = circular(0),
     to = circular(2 * pi), n = 512, K = NULL, min.k=10, control.circular=list(), \dots)
\method{print}{density.circular}(x, digits = NULL, \dots)
}

\arguments{
  \item{x}{the data from which the estimate is to be computed. The object is coerced to class \code{circular}.}
  \item{z}{the points where the density is estimated. If \code{NULL}
    equally spaced points are used according to the parameters
    \code{from}, \code{to} and \code{n}.}
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
  \item{from, to}{the left and right-most
    points of the grid at which the density is to be estimated. The objects are coerced to class \code{circular}.}
  \item{n}{the number of equally spaced points at which the density is
          to be estimated.}
  \item{K}{number of terms to be used in approximating the density.}
  \item{min.k}{minimum number of terms used in approximating the density.}
  \item{control.circular}{the attribute of the resulting objects (\code{x} component).}
  \item{digits}{integer indicating the precision to be used.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\value{
  an object with class \code{"density.circular"} whose
  underlying structure is a list containing the following components.
  \item{data}{original dataset.}
  \item{x}{the \code{n} coordinates of the points where the density is
    estimated. It is a circular objects with coordinate system setting using \code{control.circular}.}
  \item{y}{the estimated density values.}
  \item{bw}{the bandwidth used.}
  \item{N}{the sample size after elimination of missing values.}
  \item{call}{the call which produced the result.}
  \item{data.name}{the deparsed name of the \code{x} argument.}
  \item{has.na}{logical, for compatibility (always FALSE).}
}

\references{
Z.D. Bai and C.R. Rao and L.C. Zhao (1988). Kernel Estimators of Density
Function of Directional Data, Journal of Multivariate Analysis, 27, 24-39.

J. Klemel\"a (2000). Estimation of densities and derivatives of
densities with directioinal data, Journal of Multivariate Analysis, 73, 18-40.

V.R. Prayag and A.P. Gore (1990). Density Estimation for Randomly
Distributed Circular Objects, Metrika, 1990, 37, 63-69.

P. Hall and G.S. Watson and J. Cabrera (1987). Kernel Density Estimation
with Spherical Data, Biometrika, 74, 4, 751--762.
}

\author{Claudio Agostinelli}

\seealso{\code{\link{plot.density.circular}} and \code{\link{lines.density.circular}}}

\examples{
x <- rvonmises(n=100, mu=circular(pi), kappa=2)
res25 <- density(x, bw=25, control.circular=list(units="degrees"))
circularp(res25$x)
plot(res25, points.plot=TRUE, xlim=c(-1.6,1))
res50 <- density(x, bw=25, adjust=2)
lines(res50, col=2)
}

\keyword{distribution}
\keyword{smooth}
