\name{unique.circular}
\alias{unique.circular}
\title{Extract Unique Elements from a circular vector}
\description{
  \code{unique.circular} returns a circular vector but with duplicate  elements removed.
}
\usage{
unique.circular(x, incomparables = FALSE, \dots)
}

\arguments{
  \item{x}{a vector or a data frame or an array or \code{NULL}.}
  \item{incomparables}{a vector of values that cannot be compared.
    Currently, \code{FALSE} is the only possible value, meaning that all
    values can be compared.}
 \item{\dots}{}
}
\details{
  This is a method for \code{circular} object.
}

\value{
  An object of the same type of \code{x}. but if
  an element is equal to one with a smaller index, it is removed.
}

\seealso{
  \code{\link{unique}} 
}

\examples{

x <- rvonmises(10, circular(0), 10)
unique(x)

}

\keyword{manip}
\keyword{logic}

