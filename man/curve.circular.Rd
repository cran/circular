\name{curve.circular}
\alias{curve.circular}
\alias{plot.function.circular}

\title{Draw Function Plots in a Circle}
\description{
  Draws a curve corresponding to the given function or expression (in
  \code{x}) over the interval \code{[from,to]} in a circle. Mainly used to 
  plot circular density functions.
}
\usage{
curve.circular(expr, from=NULL, to=NULL, n=101, add=FALSE, 
  cex=1, axes=TRUE, ticks=FALSE, shrink=1, tcl=0.025, 
  tcl.text=0.125, tol=0.04, uin=NULL, xlim=c(-1, 1), 
  ylim=c(-1, 1), digits=2, modulo=c("2pi", "asis", "pi"), 
  main=NULL, sub=NULL, xlab="", ylab="", 
  control.circle=circle.control(), \dots)
plot.function.circular(x, from=0, to=2*pi, \dots)
}
\arguments{
  \item{expr}{an expression written as a function of \code{x}, or
    alternatively the name of a function which will be plotted.}
  \item{x}{a \sQuote{vectorizing} numeric \R function.}
  \item{from,to}{the range over which the function will be plotted.}
  \item{n}{integer; the number of x values at which to evaluate.}
  \item{add}{logical; if \code{TRUE} add to already existing plot.}
  \item{axes}{logical: if \code{TRUE} axis are added to the plot.}
  \item{ticks}{logical: if \code{TRUE} tick - marks are added to the plot.}
  \item{shrink}{parameter that controls the size of the plotted circle.
   Default is 1.  Larger values shrink the circle, while smaller values
   enlarge the circle.}
  \item{tcl}{length of the ticks.}
  \item{tcl.text}{The position of the axis labels.}
  \item{tol}{proportion of white space at the margins of plot.}
  \item{uin}{desired values for the units per inch parameter. If of length
  1, the desired units per inch on the x axis.}
  \item{xlim, ylim}{the ranges to be encompassed by the x and y
   axes. Useful for centering the plot.}
  \item{digits}{number of digits used to print axis values.}
  \item{modulo}{the modulo used to process the data.}
  \item{main, sub, xlab, ylab, cex}{graphical parameters.}
  \item{control.circle}{parameters passed to \code{\link{plot.default}} in order to draw the circle. The function \code{\link{circle.control}} is used to set the parameters.}
  \item{\dots}{parameters, passed to \code{\link{lines.circular}}.}
}
\details{
For now, curve circular draws functions defined in radians, counterclockwise coordinate and zero at 0.
}

\value{
A list with information on the plot: zero, rotation and next.points.
}

\seealso{
\code{\link{lines.circular}} and \code{\link{circle.control}}
}

\author{Claudio Agostinelli}

\examples{

ff <- function(x) sqrt(x)/20
curve.circular(ff)
curve.circular(ff, to=6*pi, join=FALSE, nosort=TRUE, n=1001, modulo="asis",
  shrink=1.2)

plot.function.circular(function(x) dvonmises(x, circular(0), 10), xlim=c(-1, 2.2))
}

\keyword{hplot}
