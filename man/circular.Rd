\name{Circular}
\alias{circular}
\alias{as.circular}
\alias{is.circular}
\alias{print.circular}
\title{Create Objects of class circular for Circular data.}
\description{
The function \code{circular} is used to create circular objects. \code{as.circular} and \code{is.circular} coerce an object to a circular and
  test whether an object is a circular data.
}
\usage{
circular(x, type = c("angles", "directions"), units = c("radians", "degrees"), template = c("none", "geographics"), modulo = c("asis", "2pi", "pi"), zero = 0, rotation = c("counter", "clock"), names)
as.circular(x, ...)
is.circular(x)
\method{print}{circular}(x, info=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a vector or a matrix.}
  \item{type}{the type of measures (Not Used Yet).}
  \item{units}{units of the measures.}
  \item{template}{how the data should be plotted (Until now only
    'geographics' style is available).}
  \item{modulo}{if we need to reduce the measures to modulo.}
  \item{zero}{the zero of the axes (in radians).}
  \item{rotation}{the orientation of the axes.}
  \item{names}{names of the data.}
  \item{info}{if \code{TRUE} information on the properties of the data
    are printed.}
  \item{...}{parameters passed to \code{circular} for \code{as.circular}
    and to \code{print.default} for \code{print.circular}.}
}

\value{
an object of class \code{\link{circular}}
}

\author{Claudio Agostinelli}

\examples{
x <- circular(c(pi, pi/3, pi/4))
print(x)
is.circular(x)
}

\keyword{misc}
