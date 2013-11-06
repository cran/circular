% File created by Eduardo Garcia-Portugues
% Modified by Claudio Agostinelli (2011/06/24)
% Modified by Claudio Agostinelli (2011/07/27)
\name{bandwidth}
\alias{bw.cv.mse.circular}
\alias{bw.cv.ml.circular}
\alias{bw.nrd.circular}
\concept{bandwidth}
\concept{Circular Bandwidth}
\title{Bandwidth Selectors for Kernel Density Estimation for Circular Data}

\description{
  Bandwidth selectors for circular kernels in \code{\link{density.circular}}.
}

\usage{
bw.cv.mse.circular(x, lower=NULL, upper=NULL, tol = 1e-4,
  kernel = c("vonmises", "wrappednormal"), K = NULL, min.k = 10)

bw.cv.ml.circular(x, lower=NULL, upper=NULL, tol = 1e-4,
  kernel = c("vonmises", "wrappednormal"), K = NULL, min.k = 10)

bw.nrd.circular(x, lower=NULL, upper=NULL,
  kappa.est=c("ML","trigmoments"), kappa.bias=FALSE, P=3)
}

\arguments{
  \item{x}{the data from which the bandwidth is to be computed. The object is coerced to class \code{circular}.}
  \item{lower, upper}{range over which to minimize for cross validatory
    bandwidths. The default is almost always satisfactory, although it
    is recommended experiment a little with different ranges. A warning
    message indicates if the resulting bandwidth is too near to the
    endpoints of the interval search.}
  \item{tol}{for cross validatory bandwidths, the convergence tolerance for \code{\link{optimize}}.}
  \item{kernel}{a character string giving the smoothing kernel to be used. This must be one of \code{"vonmises"} or \code{"wrappednormal"}.}
  \item{K}{number of terms to be used in approximating the wrappednormal density. See \code{\link{dwrappednormal}}.}
  \item{min.k}{minimum number of terms used in approximating the
    wrappednormal density. See \code{\link{dwrappednormal}}.}
  \item{kappa.est}{a numerical value or one available method.}
  \item{kappa.bias}{logical. If \code{TRUE}, when \code{kappa.est=="ML"}
    a bias correction in the estimation of kappa is applied.}
  \item{P}{integer, the maximum order of the sample trigonometric
    moments used in the estimation of \code{kappa} when
    \code{kappa.est=="trigmoments"}, see Details.}
}

\details{

\code{bw.cv.mse.circular} and \code{bw.cv.ml.circular} implement cross validatory bandwidths minimizing squared--error loss and Kullback--Leibler loss, respectively. This is done by minimizing the second and third equations in section 5 of Hall, Watson and Cabrera (1987). Kullback--Leibler loss is equivalent to maximize the cross validation log--likelihood with respect to the bandwidth parameter.

\code{bw.nrd.circular} implements a rule-of-thumb for choosing the
bandwidth of a von Mises kernel density estimator with underlaying
population von Mises. It was proposed by Taylor (2008, equation (7)) and
is the circular analogue of the usual rule of thumb used for the normal
distribution. The only remarkable difference between them is that
Taylor's bandwidth supposes a von Mises population for the derivation of
AMISE, while normal rule of thumb only introduces distribution
assumption to compute the density curvature. Estimation of the spread is
done by maximum likelihood. The "trigmoments" method for the estimation of
\code{kappa} is implemented as follows. Let \eqn{\mu_p} be the p-th
sample trigonometric moment. Let \eqn{k_p} be the estimates of
\code{kappa} using the p-th sample trigonometric moment, as solution
(using \code{uniroot} function) of the equation \eqn{A_p(k) = \frac{1}{n}
  \sum_{i=1}^n \cos(p x_i - \mu_p)}. We let \code{kappa} equal to
\eqn{max(k_1, k_2, \cdots, k_P)}, see Taylor (2008) for further details.

Note that circular bandwidth has a differet scale from linear bandwidth (see Hall, Watson and Cabrera (1987)). The behaviour of the circular bandwidth is the \emph{inverse} of the linear: large values overestimate the density, whereas small values underestimate. 
}

\section{Warning}{
Plug-in bandwidth selector \code{bw.nrd.circular} assumes that the
underlaying population is von Mises. If this is not true, it might lead
to serious misestimations of the circular bandwidth. Example 2 below
shows how this behaviour can appear with multimodality populations. In
those cases, the use of \code{kappa.est="trigmoments"} could be of help.
}

\value{
  A bandwidth on a scale suitable for the \code{bw} argument
  of \code{density.circular}.
}

\references{
P. Hall and G.S. Watson and J. Cabrera (1987). Kernel Density Estimation
with Spherical Data, Biometrika, 74, 4, 751--762.

C.C Taylor (2008). Automatic bandwidth selection for circular density estimation. Computational Statistics and Data Analysis, 53, 7, 3493--256.
}

\author{Claudio Agostinelli and Eduardo Garcia--Portugues}

\seealso{\code{\link{density.circular}}}

\examples{
set.seed(12345)

## Example 1: von Mises ##
theta1 <- rvonmises(n=150,mu=circular(pi),kappa=2)

bw.nrd1 <- bw.nrd.circular(theta1)
bw.cv.mse1 <- bw.cv.mse.circular(theta1)
bw.cv.ml1 <- bw.cv.ml.circular(theta1)

## Linear plot
plot(function(x) dvonmises(circular(x), mu=circular(pi), kappa=2),
type="l", lwd=2, col=1, main="von Mises", xlab=expression(theta),
ylab="Density", from=0, to=2*pi)
plot(approxfun(density.circular(x=theta1, bw=bw.nrd1)), col=2, from=0, to=2*pi, add=TRUE)
plot(approxfun(density.circular(x=theta1, bw=bw.cv.mse1)), col=3,
from=0, to=2*pi, add=TRUE)
plot(approxfun(density.circular(x=theta1, bw=bw.cv.ml1)), col=4, from=0,
to=2*pi, add=TRUE)
legend("topright", legend=c("True", "Taylor", "LSCV", "MLCV"), col=1:4, lwd=2)
rug(theta1)

## Circular plot
dvonmises1 <- function(x) dvonmises(circular(x), mu=circular(pi), kappa=2)
curve.circular(dvonmises1, lwd=2, col=1, main="von Mises", xlim=c(-1.5,
1.5), ylim=c(-1.5,1.5))
lines(density.circular(x=theta1, bw=bw.nrd1), col=2)
lines(density.circular(x=theta1, bw=bw.cv.mse1), col=3)
lines(density.circular(x=theta1, bw=bw.cv.ml1), col=4)
legend("topright", legend=c("True", "Taylor", "LSCV", "MLCV"), col=1:4, lwd=2)
points(theta1)

## Example 2: mixture of von Mises ##

theta2 <- rmixedvonmises(n=150, mu1=circular(pi/2),
mu2=circular(3*pi/2), kappa1=5, kappa2=5,p=0.5)

bw.nrd2 <- bw.nrd.circular(theta2)
bw.cv.mse2 <- bw.cv.mse.circular(theta2)
bw.cv.ml2 <- bw.cv.ml.circular(theta2)

## Linear plot
plot(function(x) dmixedvonmises(circular(x), mu1=circular(pi/2),
mu2=circular(3*pi/2), kappa1=5, kappa2=5, p=0.5), type="l", lwd=2,
col=1, main="mixture of von Mises", xlab=expression(theta),
ylab="Density", from=0, to=2*pi)
lines(density.circular(x=theta2, bw=bw.nrd2), plot.type='line', col=2)
lines(density.circular(x=theta2, bw=bw.cv.mse2), plot.type='line',
col=3)
lines(density.circular(x=theta2, bw=bw.cv.ml2), plot.type='line', col=4)
rug(theta2)
legend("topright", legend=c("True", "Taylor", "LSCV", "MLCV"), col=1:4, lwd=2)

## Circular plot
dmixedvonmises1 <- function(x) dmixedvonmises(circular(x), mu1=circular(pi/2),
mu2=circular(3*pi/2), kappa1=5, kappa2=5, p=0.5)
curve.circular(dmixedvonmises1, join=TRUE,
xlim=c(-1.5, 1.5), ylim=c(-1.5, 1.5), lwd=2, col=1, main="mixture of von
Mises")
lines(density.circular(x=theta2, bw=bw.nrd2), col=2)
lines(density.circular(x=theta2, bw=bw.cv.mse2), col=3)
lines(density.circular(x=theta2, bw=bw.cv.ml2), col=4)
points(theta2)
legend("topright", legend=c("True", "Taylor", "LSCV", "MLCV"), col=1:4, lwd=2)

## Example 3: mixture of von Mises and Wrapped Cauchy ##

rmixture <- function(n){
  x <- circular(sapply(runif(n), function(u) ifelse(u>0.5,
  rvonmises(n=1, mu=circular(pi),kappa=10),
  rwrappedcauchy(n=1,mu=circular(pi/2),rho=0.75))))
  return(x)
}

theta3 <- rmixture(n=150)

bw.nrd3 <- bw.nrd.circular(theta3)
bw.cv.mse3 <- bw.cv.mse.circular(theta3, lower=0.1, upper=100)
bw.cv.ml3 <- bw.cv.ml.circular(theta3, lower=0.1, upper=100)

dmixture <- function(x) (dvonmises(x, mu=circular(pi),
kappa=10)+dwrappedcauchy(x, mu=circular(pi/2), rho=0.75))/2
curve.circular(dmixture, join=TRUE, xlim=c(-1.5, 1.5), ylim=c(-1.5,
1.5), lwd=2, col=1, main="mixture of von Mises and Wrapped Normal")
lines(density.circular(x=theta3, bw=bw.nrd3), col=2)
lines(density.circular(x=theta3, bw=bw.cv.mse3), col=3)
lines(density.circular(x=theta3, bw=bw.cv.ml3), col=4)
legend("topright", legend=c("True", "Taylor", "LSCV", "MLCV"), col=1:4, lwd=2)
points(theta3)
}

\keyword{distribution}
\keyword{smooth}