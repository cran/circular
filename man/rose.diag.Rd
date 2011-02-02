\name{rose.diag}
\title{Rose Diagram}
\alias{rose.diag}

\description{
Creates a rose diagram of a circular data set on the current graphics device.}

\usage{rose.diag(x, pch = 16, cex = 1, axes = TRUE, shrink = 1, 
  bins = NULL, upper = TRUE, ticks = TRUE, tcl = 0.025, tcl.text = 0.125,
  radii.scale = c("sqrt", "linear"), border=NULL, col=NULL, tol = 0.04,
  uin = NULL, xlim = c(-1, 1), ylim = c(-1, 1), prop = 1, digits = 2, 
  plot.info = NULL, units = NULL, template = NULL, zero = NULL, 
  rotation = NULL, main = "", xlab = "", ylab = "", add = FALSE, ...)
}

\arguments{
 \item{x}{a vector, matrix or data.frame. The object is coerced to class \code{\link{circular}}.}
 \item{pch}{point character to use.  See help on \code{\link{par}}.}
 \item{cex}{point character size.  See help on \code{\link{par}}.}
 \item{axes}{logical: if \code{TRUE} axes are plotted according to properties of \code{x}.}
 \item{shrink}{parameter that controls the size of the plotted circle.
   Default is 1.  Larger values shrink the circle, while smaller values enlarge the circle.}
 \item{bins}{number of arcs to partition the circle with.}
 \item{upper}{therose diagram cells are "upper"-closed intervals.}
 \item{ticks}{logical: if \code{TRUE} ticks are plotted according to the
  value of \code{bins}.}
 \item{tcl}{length of the ticks.}
 \item{tcl.text}{the position of the axis labels.}
 \item{radii.scale}{make possible to choose sector radius form:
 square-root of relative frequency (\code{sqrt}, default) or
 conventional linear scale (\code{linear}).}
 \item{border}{the color to draw the border. The default, \code{NULL}, means to
    use \code{par("fg")}.  Use \code{border = NA} to omit borders.}
 \item{col}{the color for filling the rose diagram. The default,
   \code{NULL}, is to leave rose diagram unfilled.
color of the points. The values are recycled if needed.}
 \item{tol}{proportion of white space at the margins of plot.}
 \item{uin}{desired values for the units per inch parameter. If of length
  1, the desired units per inch on the x axis.}
 \item{xlim, ylim}{the ranges to be encompassed by the x and y
  axes. Useful for centering the plot.}
 \item{prop}{numerical constant determining the radii of the sectors.  By default, prop = 1.}
 \item{digits}{number of digits used to print axis values.}
  \item{plot.info}{an object from \code{\link{plot.circular}} that
  contains information on the \code{zero}, the \code{rotation} and \code{next.points}.}
  \item{units}{the \code{units} used in the plot. If \code{NULL} the \code{units} of the first component of 'x' is used.}
  \item{template}{the template of the plot. Ignored if \code{plot.info} is provided.}
  \item{zero}{the zero of the plot. Ignored if \code{plot.info} or \code{template} are provided.}
  \item{rotation}{the rotation of the plot. Ignored if \code{plot.info} or \code{template} are provided.}
 \item{main, xlab, ylab}{title, x label and y label of the plot.}
 \item{add}{add the rose diag to an existing plot.}
 \item{...}{futher parameters passed to \code{\link{polygon}}.}
}

\value{
a list with information on the plot: zero, rotation and next.points.
}


\note{some codes from \code{\link{eqscplot}} in MASS is used. Since
  version 0.4-1 the meaning of the \code{col} parameter is changed.}

\details{
The circumference of the circle is split into groups, the number of groups specified by bins. For each group, a sector is drawn. The radii of the sectors are by default equal to the square root of the relative frequencies of observations in each group. This ensures that the area of the sector is proportional to the group frequency. The length of the radii can be controlled by varying the parameter prop. Since version 0.3-9 the intervals are on the form [a,b).
}

\author{Claudio Agostinelli, Ulric Lund and Hiroyoshi Arai}

\seealso{
\code{\link{plot.circular}}
}

\examples{

# Generate uniform data and create several rose diagrams.  
# Some optional parameters may be needed to optimize plots.
x <- circular(runif(50, 0, 2*pi))
rose.diag(x, bins = 18, main = 'Uniform Data')
points(x)

# Generate von Mises data and create several rose diagrams.
x <- rvonmises(n=50, mu=circular(0), kappa=5, control.circular=list(zero=pi/4))
y <- rose.diag(x, bins=18) # Points fall out of bounds.
points(x, plot.info=y, stack=TRUE)
y <- rose.diag(x, bins=18, prop=1.5, shrink=1.5) # Adjust optional parameters to fit
######## all points on plot.
points(x, plot.info=y, stack=TRUE)

# Add the rose diag to a plot
plot(x)
rose.diag(x, bins=12, add=TRUE, col=2)
}



\keyword{hplot}
