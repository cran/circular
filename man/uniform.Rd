\name{Circular Uniform}
\title{Circular Uniform Density Function}
\alias{dcircularuniform}
\alias{rcircularuniform}

\description{
  Density and random generation for the Circular Uniform distribution on the whole circle. 
}

\usage{
dcircularuniform(x)
rcircularuniform(n, control.circular=list())
}

\arguments{
  \item{x}{a vector. The object is not coerced to class \code{\link{circular}}.}
\item{n}{number of observations.}
\item{control.circular}{the attribute of the resulting object.}
}

\value{
\code{dcircularuniform} gives the density and \code{rcircularuniform} generates random deviates.
}

\author{Claudio Agostinelli}

\references{
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 2.2.1, World Scientific Press, Singapore.
}

\examples{

data1 <- rcircularuniform(100, control.circular=list(units="degrees"))
plot(data1)

curve.circular(dcircularuniform, join=TRUE, xlim=c(-1.2, 1.2), 
  ylim=c(-1.2, 1.2), main="Density of a Circular Uniform Distribution")
}


\keyword{distribution}
