\name{circularp}
\alias{circularp}
\alias{circularp<-}

\title{Attributes for a Circular Object}
\description{
`circularp' returns the `circularp' attribute (or `NULL').  `circularp<-' sets the `circularp' attribute.
}
\usage{
circularp(x)
circularp(x) <- value
}

\arguments{
  \item{x}{a vector or a matrix of circular data.}
  \item{value}{a vector of length 6.}
}

\details{
     The \code{circularp} attribute is a list of five elements: type, units,
     template, modulo, zero and rotation; see \code{\link{circular}} for
     their meaning.

     Assignments are checked for consistency.    

     Assigning \code{NULL} removes the \code{circularp} attribute
     and any \code{"circular"} class of \code{x}.
}
\author{Claudio Agostinelli}

\seealso{\code{\link{circular}}}

\keyword{misc }
