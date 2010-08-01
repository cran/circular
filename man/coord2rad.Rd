\name{coord2rad}
\alias{coord2rad}
\title{
  Angles between a vector and the x-axis
}
\description{
  From coordinates of the end point of a vector in 2 dimensions to the angle between this vector and the x-axis
}
\usage{
coord2rad(x, y = NULL, control.circular = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a \code{matrix} or a \code{data.frame} with two columns if y
    is \code{NULL} otherwise a vector.}
  \item{y}{a vector.}
  \item{control.circular}{the attribute of the resulting object.}
}
\value{
an object of class \code{\link{circular}}  
}
\author{
Claudio Agostinelli and Frederick T. Wehrle  
}

\seealso{
\code{\link{circular}}  
}
\examples{
set.seed(1234)
x <- cbind(rnorm(20), rnorm(20))
y <- coord2rad(x)
}
