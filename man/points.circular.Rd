\name{points.circular}
\title{Add Points to a Circular Plot}
\alias{points.circular}

\description{
Add points to a plot of circular data points on the current graphics device. 
}

\usage{
points.circular(x, pch = 16, cex = 1, stack = FALSE, sep = 0.025, shrink = 1, 
    bins, col, next.points, plot.info, zero, rotation, ...)
}

\arguments{
\item{x}{a vector, matrix or data.frame. The object is coerced to class \code{\link{circular}}.}
\item{pch}{point character to use.  See help on \code{\link{par}}.}
\item{cex}{point character size.  See help on par.}
\item{stack}{logical: if \code{TRUE}, points are stacked on the perimeter of the circle.  Otherwise, all points are plotted on the perimeter of the circle.  Default is \code{FALSE}.}
\item{sep}{constant used to specify the distance between stacked points,
  if \code{stack==TRUE} or in the case of more than one dataset.  Default is 0.025; smaller values will create smaller spaces.}
\item{shrink}{parameter that controls the size of the plotted circle.
  Default is 1.  Larger values shrink the circle, while smaller values
  enlarge the circle.}
\item{bins}{if \code{stack==TRUE}, bins is the number of arcs to partition the circle with.}
\item{col}{color of the points. The values are recycled if needed.}
\item{next.points}{if \code{stack=FALSE}, the distance from the circle
  the next dataset is plotted. This overrides the information in \code{plot.info}.}
\item{plot.info}{an object from \code{\link{plot.circular}} that
  contains information on the \code{zero}, the \code{rotation} and \code{next.points}.}
\item{zero}{the zero of the plot. This overrides the information in \code{plot.info}.}
\item{rotation}{the rotation of the plot. This overrides the information in \code{plot.info}.}
\item{...}{futher parameters passed to \code{\link{plot.default}}.}
}

\details{
When there are many closely distributed observations, stacking is
recommended. When stacking the points, if there are many points in a particular bin, it may be necessary to shrink the plot of the circle so that all points fit.  This is controlled with the parameter \code{shrink}.  Generally the parameter \code{sep} does not need adjustment, however, when shrinking the plot, or for a very large number of observations, it may be helpful.
}

\author{Claudio Agostinelli}

\seealso{
\code{\link{plot.circular}}.
}

\examples{
data.1 <- rvonmises(n=100, mu=0, kappa=3)
data.2 <- rvonmises(n=100, mu=pi/3, kappa=3) 
res <- plot(data.1, stack=FALSE, bins=150, shrink=1.5, col=1) 
points(data.2, plot.info=res, col=2)
}

\keyword{hplot}
