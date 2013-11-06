\name{fisherB20}
\alias{fisherB20}
\alias{fisherB20c}

\title{B.20 Movements of blue periwinkles.}

\usage{
data(fisherB20)
data(fisherB20c)
}

\description{
Distances 'x' and directions 'theta' by small blue periwinkles, Nodilittorina unifasciata, after they had been transplanted downshore from the height at which they normaly live.
}
\format{
\code{fisherB20} is a \code{\link{data.frame}} of integer value.
\code{fisherB20c} is a \code{\link{data.frame}} that contains the same observations, but in the first column, the data is a \code{\link{circular}} object.
}

\source{
  N.I. Fisher (1993) pag. 252-253. Data kindly supplied by Dr A. Underwood and Ms G. Chapman.
}
\references{
  N.I. Fisher (1993) Statistical analysis of circular data. Cambridge University Press.
}

\examples{
data(fisherB20)
data(fisherB20c)
par(mfcol=c(1,3))
plot(fisherB20c$theta, xlab=expression(theta))
boxplot(fisherB20c$x, xlab="x")
plot(c(fisherB20$x, fisherB20$x), c(fisherB20$theta,
  fisherB20$theta+360), xlab="x", ylab=expression(theta))

}

\keyword{datasets}
