\name{lsfit.circle}
\alias{lsfit.circle}
\alias{print.lsfit.circle}

\title{Fit a 2D circle to an (x,y) dataset}
\description{
Fit a 2D circle to an (x,y) dataset using LS.
}
\usage{
lsfit.circle(x, y, init = NULL, units = c("radians", "degrees"), 
  template = c("none", "geographics"),
  modulo = c("asis", "2pi", "pi"), zero = 0, 
  rotation = c("counter", "clock"), ...)
\method{print}{lsfit.circle}(x, digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{x}{either a matrix with two columns or a vector.}
  \item{y}{if \code{x} is a vector then \code{y} must be a vector with the same length.}
  \item{init}{initial values of the parameters. A vector of length 3 with
  the following components: radius of the circle, x-coordinate of the
  center, y-coordinate of the center. If \code{NULL} the vector is set
  to \code{c(max(c(abs(x-mean(x)), abs(y-mean(y)))), mean(x), mean(y)}.}
  \item{units}{the \code{units} used in defining the angles between
    observations and the center of the circle. See \code{\link{circular}}.}
  \item{template}{the \code{template} used in defining the angles between
    observations and the center of the circle. See \code{\link{circular}}.}
  \item{modulo}{the \code{modulo} used in defining the angles between
    observations and the center of the circle. See \code{\link{circular}}.}
  \item{zero}{the \code{zero} used in defining the angles between
    observations and the center of the circle. See \code{\link{circular}}.}
  \item{rotation}{the \code{rotation} used in defining the angles between
    observations and the center of the circle. See \code{\link{circular}}.}
  \item{\dots}{further parameters passed to the \code{optim} function.}
  \item{digits}{the number of digits to be printed.}  
}

\details{
\code{lsfit.circle} uses the \code{optim} function to minimize the sum of the
squared residuals between the observations and the optimally fitting circle.
}

\value{

  An object of class \code{lsfit.circle}.

  \item{coefficients}{a vector of length 3 with the estimated radius and
  coordinate of the center of the fitted circle.}
  \item{x}{the x-coordinate.}
  \item{y}{the y-coordinate.}
  \item{x.centered}{the x-coordinate re-centered at the center of the circle.}
  \item{y.centered}{the y-coordinate re-centered at the center of the circle.}
  \item{angles}{angles of the observations with respect to the center
    coordinate of class \code{circular}.}
  \item{radius}{the distance between the observations and the center
    coordinate}
  \item{convergence}{value from the function optim.}
  \item{optim}{the output from the function optim.}
  \item{call}{the \code{\link[base]{match.call}} result.}
}
\references{ 
Coope, I. (1993). Circle fitting by linear and non-linear least squares. 
Journal of Optimization Theory and Applications, 76, 381-388. 
}
\author{Claudio Agostinelli and Ulric Lund}

\examples{

data(coope)
res <- lsfit.circle(x=x.coope, y=y.coope)
res

plot(res)

par(mfcol=c(1,2))
plot(res$angles)
hist(res$radius)

plot(circular(0), type="n", xlim=c(-5.2, 5.2), ylim=c(-5.2, 5.2), 
  xlab="The Radius of the circle \n is measured from the base line of the axes.")
lines(x=res$angles, y=res$radius, join=TRUE, type="b")
ff <- function(x) sqrt((res$coefficients[1]*cos(x))^2+(res$coefficients[1]*sin(x))^2)
curve.circular(ff, add=TRUE, join=TRUE, nosort=FALSE, col=2)

windrose(x=res$angles, y=res$radius)
}

\keyword{models}
