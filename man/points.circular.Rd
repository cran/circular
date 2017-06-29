\name{points.circular}
\title{Add Points to a Circular Plot}
\alias{points.circular}

\description{
Add points to a plot of circular data points on the current graphics device. 
}

\usage{
\method{points}{circular}(x, pch = 16, cex = 1, stack = FALSE,
  start.sep=0, sep = 0.025, 
  shrink = 1, bins = NULL, col = NULL, next.points = NULL, 
  plot.info = NULL, zero = NULL, rotation = NULL, ...)
}

\arguments{
\item{x}{a vector, matrix or data.frame. The object is coerced to class \code{\link{circular}}.}
\item{pch}{point character to use.  See help on \code{\link{par}}.}
\item{cex}{point character size.  See help on par.}
\item{stack}{logical: if \code{TRUE}, points are stacked on the
  perimeter of the circle.  Otherwise, all points are plotted on the
  perimeter of the circle.  Default is \code{FALSE}.}
\item{start.sep}{constant used to specify the distance between the
  center of the point and the axis.}
\item{sep}{constant used to specify the distance between stacked points,
  if \code{stack==TRUE} or in the case of more than one dataset.  Default is 0.025; smaller values will create smaller spaces.}
\item{shrink}{parameter that controls the size of the plotted circle.
  Default is 1.  Larger values shrink the circle, while smaller values
  enlarge the circle.}
\item{bins}{if \code{stack==TRUE}, bins is the number of arcs to partition the circle with.}
\item{col}{color of the points. The values are recycled if needed.}
\item{next.points}{if \code{stack=FALSE}, the distance from the circle
  the next dataset is plotted. Ingored if \code{plot.info} is provided.}
\item{plot.info}{an object from \code{\link{plot.circular}} that
  contains information on the \code{zero}, the \code{rotation} and \code{next.points}.}
\item{zero}{the zero of the plot. Ignored if \code{plot.info} is provided.}
\item{rotation}{the rotation of the plot. Ignored if \code{plot.info} is provided.}
\item{...}{futher parameters passed to \code{\link{points.default}}.}
}

\details{
When there are many closely distributed observations, stacking is
recommended. When stacking the points, if there are many points in a particular bin, it may be necessary to shrink the plot of the circle so that all points fit. This is controlled with the parameter \code{shrink}. Generally the parameter \code{sep} does not need adjustment, however, when shrinking the plot, or for a very large number of observations, it may be helpful. Since version 0.3-9 the intervals are on the form [a,b).
}

\author{Claudio Agostinelli}

\value{
A list with information on the plot: zero, rotation and next.points.
}

\seealso{
\code{\link{plot.circular}} and \code{\link{lines.circular}}.
}

\examples{
data.1 <- rvonmises(n=100, mu=circular(0), kappa=3)
data.2 <- rvonmises(n=100, mu=circular(pi/3), kappa=3) 
res <- plot(data.1, stack=FALSE, col=1) 
points(data.2, plot.info=res, col=2)
}

\keyword{hplot}
