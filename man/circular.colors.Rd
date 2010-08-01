\name{circular.colors}
\alias{circular.colors}

\title{
Color Palettes for Circular
}

\description{
Create a vector of \code{n} contiguous colors.
}

\usage{
circular.colors(n, m = 0, M = 2 * pi, offset = 0, ...)
}

\arguments{
  \item{n}{the number of colors (>= 1) to be in the palette.}
  \item{m}{the smallest angle in radians.}
  \item{M}{the largest angle in radians.}
  \item{offset}{the zero in radians.}
  \item{\dots}{further arguments passed to the function \code{\link{hsv}}.}
}

\value{
  a vector of length \code{n}.
}

\author{
  Claudio Agostinelli
}

\seealso{
  \code{\link{hsv}}, \code{\link{colors}}
}

\examples{
  circular.colors(n=10, m=0, M=2*pi)
}

\keyword{color}
