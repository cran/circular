\name{c.circular}
\alias{c.circular}

\title{A method for circular object, which combines its arguments}
\description{

}

\usage{
c.circular(..., recursive = FALSE)
}

\arguments{
  \item{\dots}{vectors, the first of which of class \code{circular}.}
  \item{recursive}{logical. If 'recursive=TRUE', the function recursively
          descends through lists combining all their elements into a
          vector.
}
}

\author{Claudio Agostinelli}

\seealso{\code{\link{c}} }

\examples{

x <- rvonmises(10, circular(0), 10)
y <- rvonmises(10, circular(0), 10, control.circular=list(units="degrees"))
z <- runif(10, 0, 20) # here you do not use circular properties, but you mean it is measured in degrees
c(x, y, z) # While y is converted in radians, z is treated as it was!

}
\keyword{manip}

