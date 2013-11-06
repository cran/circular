\name{asytriangular}    
\title{Asymmetric Triangular Density Function}
\alias{dasytriangular}

\description{
Density the Asymmetric Triangular circular distribution. 
}

\usage{
dasytriangular(x, rho)
}

\arguments{
\item{x}{a vector. The object is coerced to class \code{\link{circular}}.}
\item{rho}{concentration parameter of the distribution.  rho must be
  between 0 and \eqn{1/pi}.}
}

\value{
\code{dasytriangular} gives the density. 
}

\author{Claudio Agostinelli}

\references{
Mardia (1972) Statistics for Directional Data, Wiley. Pag. 52
}

\examples{

ff <- function(x) dasytriangular(x, rho=0.3)
curve.circular(ff, shrink=1.2, join=TRUE)

}



\keyword{distribution}
