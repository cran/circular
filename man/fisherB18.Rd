\name{fisherB18}
\alias{fisherB18}
\alias{fisherB18c}

\title{B.18 Wind direction and ozone concentration.}

\usage{
data(fisherB18)
data(fisherB18c)
}

\description{
19 measuremenents of wind direction 'theta' and ozone level 'x' taken at 6.00am at four-day intervals between April 18th and June 29th, 1975 at a weather station in Milwaukee.
}
\format{
\code{fisherB18} is a \code{\link{data.frame}} of integer value.
\code{fisherB18c} is a \code{\link{data.frame}} that contains the same observations, but in the first column, the data is a \code{\link{circular}} object.
}

\source{
  N.I. Fisher (1993) pag. 251. Johnson & Wehrly (1977, Table 1).
}
\references{
  N.I. Fisher (1993) Statistical analysis of circular data. Cambridge University Press.
}

\examples{
data(fisherB18)
data(fisherB18c)
par(mfcol=c(1,3))
plot(fisherB18c$theta, xlab=expression(theta))
boxplot(fisherB18c$x, xlab="x")
plot(c(fisherB18$x, fisherB18$x), c(fisherB18$theta,
  fisherB18$theta+360), xlab="x", ylab=expression(theta))
}

\keyword{datasets}
