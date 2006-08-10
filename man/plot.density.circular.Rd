\name{plot.density.circular}
\alias{plot.density.circular}

\title{Plot Method for Kernel Density Estimation for Circular Data}
\description{
   The \code{plot} method for \code{density.circular} objects.
}
\usage{
\method{plot}{density.circular}(x, main=NULL, xlab=NULL, ylab="Density circular", type="l",
           zero.line=TRUE, points.plot=FALSE, points.col=1, points.pch=1, points.cex=1,
           plot.type=c("circle", "line"), axes=TRUE, ticks=TRUE, bins=NULL, shrink=1, 
           tcl=0.025, tcl.text = 0.125, sep=0.025, tol=0.04, digits=2, cex=1, uin=NULL, 
           xlim=NULL, ylim=NULL, join=FALSE, nosort=FALSE, units=NULL, 
           template=NULL, zero=NULL, rotation=NULL, ...)
}

\arguments{
  \item{x}{an object of class \code{\link{density.circular}}.}
  \item{main, xlab, ylab, type}{plotting parameters with useful defaults.}
  \item{zero.line}{logical; if \code{TRUE}, add a base line at \eqn{y =
      0}. Used when \code{plot.type} is \code{"line"}.}
  \item{points.plot}{logical; if \code{TRUE} original data are added to
    the plot.}
  \item{points.col, points.pch, points.cex}{parameters used to draw the points.}
  \item{plot.type}{type of the plot: "line": linear plot, "circle": circular plot.}
  \item{axes}{logical; if \code{TRUE} axis are drawn.}
  \item{ticks}{logical; if \code{TRUE} ticks are drawn.}
  \item{bins}{number of ticks to plot. }
  \item{shrink}{parameter that controls the size of the plotted circle.
  Default is 1.  Larger values shrink the circle, while smaller values
  enlarge the circle.}
 \item{tcl}{length of the ticks.}
 \item{tcl.text}{The position of the axis labels.}
 \item{sep}{constant used to specify the distance between stacked points. Default is 0.025; smaller values will create smaller spaces.}
 \item{tol}{proportion of white space at the margins of plot}
 \item{digits}{number of digits used to print axis values.}
 \item{cex}{point character size.  See help on \code{\link{par}}.}
 \item{uin}{desired values for the units per inch parameter. If of length
  1, the desired units per inch on the x axis.}
 \item{xlim, ylim}{the ranges to be encompassed by the x and y axes. Useful for centering the plot.}
 \item{join}{logical: should the first and the last point joined.}
 \item{nosort}{logical: should the data sort before plotting. Defaults is to sort.}
 \item{units}{units measure used in the plot. If \code{NULL} the value is taken from the attribute of object 'x' from the argument 'x', i.e. \code{x$x}.}
 \item{template}{template used in the plot.  If \code{NULL} the value is taken from the attribute of object 'x' from the argument 'x', i.e. \code{x$x}.}
 \item{zero}{position of the zero used in the plot.  If \code{NULL} the value is taken from the attribute of object 'x' from the argument 'x', i.e. \code{x$x}.}
 \item{rotation}{rotation used in the plot.  If \code{NULL} the value is taken from the attribute of object 'x' from the argument 'x', i.e. \code{x$x}.}
 \item{\dots}{futher parameters passed to \code{\link{plot.default}}.}
}

\value{
A list with information on the plot: zero, rotation and next.points.
}

\author{Claudio Agostinelli}

\seealso{\code{\link{density.circular}}, \code{\link{lines.density.circular}}, \code{\link{plot.circular}}, \code{\link{lines.circular}} and \code{\link{curve.circular}}.}

\examples{

set.seed(1234)
x <- rvonmises(n=100, mu=circular(pi), kappa=2)
res25x <- density(x, bw=25)
plot(res25x, points.plot=TRUE, xlim=c(-1.5,1))
res50x <- density(x, bw=25, adjust=2)
lines(res50x, col=2)

resp25x <- plot(res25x, points.plot=TRUE, xlim=c(-1, 1.3), ylim=c(-1.5,1.2), template="geographics", main="Plotting density estimate for two data set")
y <- rvonmises(n=100, mu=circular(pi/2), kappa=2, control.circular=list(template="geographics"))
res25y <- density(y, bw=25)
lines(res25y, points.plot=TRUE, plot.info=resp25x, col=2, points.col=2)

plot(res25x, plot.type="line", points.plot=TRUE, xlim=c(-1, 1.3), ylim=c(-1.5,1.2), template="geographics", main="Plotting density estimate for two data set")
lines(res25y, plot.type="line", points.plot=TRUE, col=2, points.col=2)

}

\keyword{dplot}

