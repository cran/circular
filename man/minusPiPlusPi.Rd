\name{minusPiPlusPi}
\alias{minusPiPlusPi}
\title{
   return angles in the [-pi,pi] interval. 
}
\description{
   return angles in the (-pi,pi] interval. 
}
\usage{
  minusPiPlusPi(x)
}
\arguments{
  \item{x}{an object of class \code{\link{circular}}.}
}
\value{
  a \code{\link{circular}} object with values in the interval (-pi,pi].
}
\author{
  Claudio Agostinelli and Alessandro Gagliardi
}

\examples{
  x <- circular(c(0, 90, 180, 270), units="degrees")
  minusPiPlusPi(x)
}

\keyword{univariate}
