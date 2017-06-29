\name{pp.unif.plot}    
\title{Uniform Circular Probability-Probability Plot}
\alias{pp.unif.plot}

\description{
Plots the empirical distribution of a data set against a uniform circular distribution function.
}

\usage{
pp.unif.plot(x, ref.line = TRUE, frac = NULL,  xlab = "Uniform Distribution", 
  ylab = "Empirical Distribution", col = NULL,
  col.inf = NULL, col.sup = NULL, ...)
}

\arguments{
\item{x}{a vector. The object is coerced to class \code{\link{circular}}.}
\item{ref.line}{logical, if TRUE a 45 degree reference line is added to
  the plot. Default is TRUE.}
\item{frac}{a number in the [0,1] interval or \code{NULL}.}
\item{xlab, ylab}{labels of the axis.}
\item{col}{color of the points.}
\item{col.inf, col.sup}{color of the \code{frac}tion of the points replicated in the left bottom and right upper corner of the plot.}    
\item{...}{parameters passed to the \code{\link{plot.default}} function.}
}

\author{Claudio Agostinelli}

\seealso{
\code{\link{pp.plot} for the von Mises distribution.}
}

\examples{
x <- rvonmises(n=25, mu=circular(0), kappa=3)
pp.unif.plot(x)
pp.unif.plot(x, frac=0.2)  
}

\keyword{hplot}
