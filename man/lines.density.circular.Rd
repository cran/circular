\name{lines.density.circular}
\alias{lines.density.circular}

\title{Add a Plot for Kernel Density Estimation for Circular Data}
\description{
   The \code{lines} add a plot for \code{density.circular} objects.
}
\usage{
\method{lines}{density.circular}(x, type = "l", zero.line = TRUE, points.plot = FALSE, points.col = 1, points.pch = 1, plot.type = c("circle", "line"), bins, shrink = 1, tcl = 0.025, ...)
}

\arguments{
  \item{x}{an object of class \code{\link{density.circular}}.}
  \item{type}{plotting parameter with useful default.}
  \item{zero.line}{logical; if \code{TRUE}, add a base line at \eqn{y =
      0}. Used when \code{plot.type} is \code{"line"}.}
  \item{points.plot}{logical; if \code{TRUE} original data are added to
    the plot.}
  \item{points.col, points.pch}{parameters used to draw the points.}
  \item{plot.type}{type of the plot.}
  \item{bins}{number of ticks to plot. }
  \item{shrink}{parameter that controls the size of the plotted circle.
  Default is 1.  Larger values shrink the circle, while smaller values
  enlarge the circle.}
\item{tcl}{length of the ticks.}
  \item{\dots}{futher parameters passed to \code{\link{lines.default}}.}
}

\value{
None.
}

\author{Claudio Agostinelli}

\seealso{\code{\link{density.circular}} and \code{\link{plot.density.circular}} }

\examples{
x <- rvonmises(n=100, mu=pi, kappa=2)
res25 <- density(x, bw=25)
plot(res25, points.plot=TRUE, xlim=c(-1.5,1))
res50 <- density(x, bw=25, adjust=2)
lines(res50, col=2)
}

\keyword{dplot}

