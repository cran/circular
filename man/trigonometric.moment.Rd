\name{trigonometric.moment}    
\title{Trigonometric Moments}
\alias{trigonometric.moment}  

\description{
Computes the specified order trigonometric moment for a set of directional data points. 
}

\usage{
trigonometric.moment(x, p=1, center=FALSE)
}

\arguments{
\item{x}{ }
\item{p}{order of trigonometric moment to be computed. Default is for a first order trigonometric moment.}
\item{center}{logical, whether to compute centered moments or not.  Default is to compute an uncentered moment.}
}

\value{
Returns a list with variables mu, rho, cos, sin, p, n, call,
respectively the pth trigonometric moment's direction, resultant length,
real and imaginary components, the order, the number of observations and
the call.
}

\author{Claudio Agostinelli and Ulric Lund}

\references{
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 1.3, World Scientific Press, Singapore.
}

\seealso{
\code{\link{var.circular}}, \code{\link{mean.circular}},
\code{\link{summary.circular}}, \code{\link{mle.vonmises}} and \code{\link{rho.circular}}
}

\keyword{math}
