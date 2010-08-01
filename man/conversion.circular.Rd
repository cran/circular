\name{conversion.circular}
\alias{conversion.circular}
\title{Unit of Measure Conversion for Circular Data and other conversions}
\description{
Conversion for Circular Data from one coordinate/units system to another one. For back compatibility, without arguments the function converts data from degrees to radians.
}

\usage{
conversion.circular(x, units = c("radians", "degrees", "hours"), type = NULL, 
  template = NULL, modulo = NULL, zero = NULL, rotation = NULL)
}

\arguments{
  \item{x}{an object of class \code{\link{circular}}.}
  \item{units}{unit of the transformed data.}
  \item{type}{type of the transformed data. If \code{NULL} no action is performed.}
  \item{template}{template of the tranformed data. If \code{NULL} no action is performed.}
  \item{modulo}{modulo of the transformed data. If \code{NULL} no action is performed.}
  \item{zero}{zero of the tranformed data. If \code{NULL} no action is performed.}
  \item{rotation}{rotation of the tranformed data. If \code{NULL} no action is performed.}
}

\value{
  an object of class \code{\link{circular}} with the specified unit of measure, modulo, zero and rotation.
}

\author{Claudio Agostinelli}

\seealso{
\code{\link{deg}} and \code{\link{rad}}. If you want to set the properties of an object instead to tranform it, you can use \code{\link{circular}} or \code{\link{circularp<-}}.
}

\examples{
x <- rvonmises(n=10, mu=circular(0), kappa=9, control.circular=list(units="degrees"))
par(mfcol=c(2, 2))
plot(x)
y <- conversion.circular(x) # only the unit is changed (to radians) and 
####### the data converted.
plot(y)
z <- conversion.circular(x, units="degrees", zero=pi) # only the zero is changed and 
####### the data converted.
plot(z)
w <- conversion.circular(x, zero=pi, rotation="clock") # zero and rotation is 
####### changed and the data converted.
plot(w)
}

\keyword{misc}
