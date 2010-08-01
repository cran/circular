\name{mixedvonmises}    
\title{Mixture of von Mises Distributions}
\alias{dmixedvonmises}
\alias{rmixedvonmises}

\description{
  Density and random generation for the mixed von Mises circular distribution.  
}

\usage{
dmixedvonmises(x, mu1, mu2, kappa1, kappa2, p)
rmixedvonmises(n, mu1, mu2, kappa1, kappa2, p, control.circular = list())
}

\arguments{
  \item{x}{a vector. The object is coerced to class \code{\link{circular}}.}
  \item{n}{number of observations.}
\item{mu1}{mean direction of one of the two von Mises distributions as a \code{circular} object.}
\item{mu2}{mean direction of the other von Mises distribution as a \code{circular} object.}
\item{kappa1}{concentration parameter of one of the two von Mises distributions.}
\item{kappa2}{concentration parameter of the other von Mises distribution.}
\item{p}{mixing proportion.}
\item{control.circular}{the attribute of the resulting object.}
}

\value{
\code{dmixedvonmises} gives the density and \code{rmixedvonmises} generates random deviates.
}

\author{Claudio Agostinelli and Ulric Lund}

\seealso{
\code{\link{dvonmises}} and \code{\link{rvonmises}}
}

\examples{

x <- rmixedvonmises(n=100, mu1=circular(0), mu2=circular(pi), kappa1=15, 
  kappa2=15, p=0.5)
plot(x)


}


\keyword{distribution}
