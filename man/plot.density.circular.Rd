\name{plot.density.circular}
\alias{plot.density.circular}

\title{Plot Method for Kernel Density Estimation for Circular Data}
\description{
   The \code{plot} method for \code{density.circular} objects.
}
\usage{
\method{plot}{density.circular}(x, main = NULL, xlab = NULL, ylab = "Density circular", type = "l", zero.line = TRUE, points.plot = FALSE, points.col = 1, points.pch = 1, plot.type = c("circle", "line"), axes = TRUE, ticks = TRUE, bins, shrink = 1, tcl = 0.025, tol = 0.04, uin, xlim = c(-1, 1), ylim = c(-1, 1), ...)
}

\arguments{
  \item{x}{an object of class \code{\link{density.circular}}.}
  \item{main, xlab, ylab, type}{plotting parameters with useful defaults.}
  \item{zero.line}{logical; if \code{TRUE}, add a base line at \eqn{y =
      0}. Used when \code{plot.type} is \code{"line"}.}
  \item{points.plot}{logical; if \code{TRUE} original data are added to
    the plot.}
  \item{points.col, points.pch}{parameters used to draw the points.}
  \item{plot.type}{type of the plot.}
  \item{axes}{logical; if \code{TRUE} axis are drawn.}
  \item{ticks}{logical; if \code{TRUE} ticks are drawn.}
  \item{bins}{number of ticks to plot. }
  \item{shrink}{parameter that controls the size of the plotted circle.
  Default is 1.  Larger values shrink the circle, while smaller values
  enlarge the circle.}
\item{tcl}{length of the ticks.}
\item{tol}{proportion of white space at the margins of plot}
\item{uin}{desired values for the units per inch parameter. If of length
  1, the desired units per inch on the x axis.}
\item{xlim, ylim}{the ranges to be encompassed by the x and y
  axes. Useful for centering the plot.}
  \item{\dots}{futher parameters passed to \code{\link{plot.default}}.}
}

\value{
None.
}

\author{Claudio Agostinelli}

\seealso{\code{\link{density.circular}} and \code{\link{lines.density.circular}} }

\examples{
x <- rvonmises(n=100, mu=pi, kappa=2)
res25 <- density(x, bw=25)
plot(res25, points.plot=TRUE, xlim=c(-1.5,1))
res50 <- density(x, bw=25, adjust=2)
lines(res50, col=2)
}

\keyword{dplot}

