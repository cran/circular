\name{rose.diag}
\title{Rose Diagram}
\alias{rose.diag}

\description{
Creates a rose diagram of a circular data set on the current graphics device.}

\usage{rose.diag(x, pch = 16, axes = TRUE, shrink = 1, bins, ticks = TRUE, tcl = 0.025, col, tol = 0.04,
uin, xlim = c(-1, 1), ylim = c(-1, 1), prop = 1, main = NULL, ...)
}

\arguments{

\item{x}{a vector, matrix or data.frame. The object is coerced to class \code{\link{circular}}.}
\item{pch}{point character to use.  See help on \code{\link{par}}.}
\item{axes}{logical: if \code{TRUE} axes are plotted according to properties of \code{x}.}
\item{shrink}{parameter that controls the size of the plotted circle.
  Default is 1.  Larger values shrink the circle, while smaller values
  enlarge the circle.}
\item{bins}{if \code{stack==TRUE}, bins is the number of arcs to partition the circle with.}
\item{ticks}{logical: if \code{TRUE} ticks are plotted according to the
  value of \code{bins}.}
\item{tcl}{length of the ticks.}
\item{col}{color of the points. The values are recycled if needed.}
\item{tol}{proportion of white space at the margins of plot.}
\item{uin}{desired values for the units per inch parameter. If of length
  1, the desired units per inch on the x axis.}
\item{xlim, ylim}{the ranges to be encompassed by the x and y
  axes. Useful for centering the plot.}
\item{prop}{numerical constant determining the radii of the sectors.  By default, prop = 1.}
\item{main}{title of the plot.}
\item{...}{futher parameters passed to \code{\link{plot.default}}.}
}

\value{
a list with information on the plot: zero, rotation and next.points.
}


\note{some codes from \code{\link{eqscplot}} in MASS is used.}

\details{
The circumference of the circle is split into groups, the number of groups specified by bins.  For each group, a sector is drawn.  The radii of the sectors are by default equal to the square root of the relative frequencies of observations in each group.  This ensures that the area of the sector is proportional to the group frequency.  The length of the radii can be controlled by varying the parameter prop. 
}

\author{Claudio Agostinelli and Ulric Lund}

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
x <- rvonmises(n=25, mu=0, kappa=5, zero=pi/4)
y <- rose.diag(x, bins=18) # Points fall out of bounds.
points(x, plot.info=y, stack=TRUE)
y <- rose.diag(x, bins=36, prop=1.5, shrink=1.5) # Adjust optional parameters to fit all points on plot.
points(x, plot.info=y, stack=TRUE)
}

\keyword{hplot}
