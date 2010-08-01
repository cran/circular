\name{lines.density.circular}
\alias{lines.density.circular}

\title{Add a Plot for Kernel Density Estimation for Circular Data}
\description{
   The \code{lines} add a plot for \code{density.circular} objects.
}
\usage{
\method{lines}{density.circular}(x, type = "l", zero.line = TRUE,
  points.plot = FALSE, points.col = 1, points.pch = 1, points.cex = 1,
  plot.type = c("circle", "line"), bins = NULL, offset=1, shrink = 1,
  tcl = 0.025, sep = 0.025, join = TRUE, nosort = FALSE,
  plot.info = NULL, zero = NULL, rotation = NULL, ...)
}

\arguments{
  \item{x}{an object of class \code{\link{density.circular}}.}
  \item{type}{plotting parameter with useful default.}
  \item{zero.line}{logical; if \code{TRUE}, add a base line at \eqn{y =
      0}. Used when \code{plot.type} is \code{"line"}.}
  \item{points.plot}{logical; if \code{TRUE} original data are added to
    the plot.}
  \item{points.col, points.pch, points.cex}{parameters used to draw the points.}
  \item{plot.type}{type of the plot.}
  \item{bins}{number of ticks to plot. }
  \item{offset}{the radius of the circle}  
  \item{shrink}{parameter that controls the size of the plotted function.
  Default is 1.}
  \item{tcl}{length of the ticks.}
  \item{sep}{constant used to specify the distance between stacked points. Default is 0.025; smaller values will create smaller spaces.}
  \item{join}{logical: should the first and the last point joined.}
  \item{nosort}{logical: should the data sort before plotting. Defaults is to sort.}
  \item{plot.info}{an object from \code{\link{plot.circular}} that
  contains information on the \code{zero}, the \code{rotation} and \code{next.points}.}
  \item{zero}{the zero of the plot. Ignored if \code{plot.info} is provided.}
  \item{rotation}{the rotation of the plot. Ignored if \code{plot.info} is provided.}
  \item{\dots}{futher parameters passed to \code{\link{lines.default}}.}
}

\value{
A list with information on the plot: zero, rotation and next.points and,
if available, the coordinates x and y.
}

\author{Claudio Agostinelli}

\seealso{\code{\link{density.circular}} and \code{\link{plot.density.circular}} }

\examples{
set.seed(1234)
x <- rvonmises(n=100, mu=circular(pi), kappa=2)
y <- rvonmises(n=100, mu=circular(pi/2), kappa=2)
resx <- density(x, bw=25)
res <- plot(resx, points.plot=TRUE, xlim=c(-1.5,1), ylim=c(-1.1, 1.5))
resy <- density(y, bw=25)
lines(resy, points.plot=TRUE, col=2, points.col=2, plot.info=res)
}

\keyword{dplot}

