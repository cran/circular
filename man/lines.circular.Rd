\name{lines.circular}
\alias{lines.circular}
\title{Add Connected Line Segments to a Circular Plot}
\description{
A method taking coordinates in a polar system and joining the corresponding points with line segments.

}
\usage{
\method{lines}{circular}(x, y, join = FALSE, nosort = FALSE, offset=1, shrink=1,
  plot.info = NULL, zero = NULL, rotation = NULL, modulo = NULL, \dots)
}

\arguments{
  \item{x}{a vector of class \code{circular}.}
  \item{y}{a vector with the same length as 'x'.}
  \item{join}{logical: if \code{TRUE} the first and the last values are joined by a line.}
  \item{nosort}{logical: if \code{TRUE} the data are not sorted before
    join them.}
  \item{offset}{the radius of the circle}
  \item{shrink}{parameter that controls the size of the plotted function.
    Default is 1.}
  \item{plot.info}{an object from another circular graphic function.}
  \item{zero}{the zero of the axis.}
  \item{rotation}{the rotation of the axis.}
  \item{modulo}{the modulo applied to 'x' before sorting.}
  \item{\dots}{graphical parameters passed to \code{\link{lines.default}}.}
}

\value{
A list with information on the plot: zero, rotation and next.points.
}

\author{Claudio Agostinelli}

\seealso{\code{\link{plot.circular}}}

\examples{

x <- rvonmises(20, circular(0), 10)
y <- runif(20, 0.5, 1)

plot(x, shrink=2)
lines(x, y)

}

\keyword{hplot}

