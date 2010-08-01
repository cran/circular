\name{circular}
\alias{circular}
\alias{as.circular}
\alias{is.circular}
\alias{print.circular}
\title{Create Objects of class circular for Circular data.}
\description{
The function \code{circular} is used to create circular objects. \code{as.circular} and \code{is.circular} coerce an object to a circular and  test whether an object is a circular data.
}
\usage{
circular(x, type = c("angles", "directions"), 
  units = c("radians", "degrees", "hours"),
  template = c("none", "geographics", "clock12", "clock24"), modulo = c("asis", "2pi", "pi"), 
  zero = 0, rotation = c("counter", "clock"), names)
as.circular(x, control.circular=list(), ...)
is.circular(x)
\method{print}{circular}(x, info=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a vector or a matrix. If a data.frame is supply then it is corced to a matrix.}
  \item{type}{the type of measures (Not Used Yet).}
  \item{units}{units of the measures.}
  \item{template}{how the data should be plotted. This set \code{modulo}, \code{zero} and \code{rotation} to some suitable values. For instance for 'geographics': zero=pi/2 and rotation='clock'. It is also used to set default labels on the plots.}
  \item{modulo}{if we need to reduce the measures to modulo.}
  \item{zero}{the zero of the axes (in radians, counter).}
  \item{rotation}{the orientation of the axes.}
  \item{names}{names of the data.}
  \item{info}{if \code{TRUE} information on the properties of the data
    are printed.}
  \item{control.circular}{the attribute (coordinate system) used to coerced the resulting objects. See \code{\link{circular}}.}
  \item{...}{For \code{as.circular} an alternative way of setting the coordinate system of the resulting objects. Passed parameters to \code{print.default} for \code{print.circular}.}
}

\value{
an object of class \code{\link{circular}}. Since version 0.3-5 the previous class of the object is retain.
}

\author{Claudio Agostinelli}

\seealso{
\code{\link{conversion.circular}}
}

\examples{
x <- circular(c(pi, pi/3, pi/4))
print(x)
is.circular(x)

x <- circular(runif(10, -pi/2, pi/2), template="geographics")
plot(x)
class(x)

x <- circular(data.frame(runif(10, -pi/2, pi/2)))
plot(x)
class(x)

cbind(x, x) # the matrix, cbind, rbind functions unclass and lost attributes! 
########Use it with care.

x <- c(pi/12,2*pi+pi/12)%%(2*pi)  # unique may not work as desidered due to machine precision 
print(x)
x <- unique(x)
print(x)

x[1]==x[2]

all.equal(x[1], x[2])

x <- as.circular(pi, control.circular=list(units="radians", zero=pi))
y <- conversion.circular(circular(pi), zero=pi)
res <- plot(x)
points(y, col=2, plot.info=res)

}

\keyword{misc}
