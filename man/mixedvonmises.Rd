\name{mixedvonmises}    
\title{Mixture of von Mises Distributions}
\alias{dmixedvonmises}
\alias{rmixedvonmises}

\description{
  Density and random generation for the mixed von Mises circular distribution.  
}

\usage{
dmixedvonmises(x, mu1, mu2, kappa1, kappa2, p)
rmixedvonmises(n, mu1, mu2, kappa1, kappa2, p, units = c("radians",
"degrees"), \dots)
}

\arguments{
  \item{x}{a vector. The object is coerced to class \code{\link{circular}}.}
  \item{n}{number of observations.}
\item{mu1}{mean direction of one of the two von Mises distributions.}
\item{mu2}{mean direction of the other von Mises distribution.}
\item{kappa1}{concentration parameter of one of the two von Mises distributions.}
\item{kappa2}{concentration parameter of the other von Mises distribution.}
\item{p}{mixing proportion.}
\item{units}{units of measure desired.}
\item{\dots}{parameters passed to \code{\link{circular}}.}
}

\value{
\code{dmixedvonmises} gives the density and \code{rmixedvonmises} generates random deviates.
}

\author{Claudio Agostinelli and Ulric Lund}

\seealso{
\code{\link{dvonmises}} and \code{\link{rvonmises}}
}

\keyword{distribution}
