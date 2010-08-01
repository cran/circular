\name{pp.plot}    
\title{von Mises Probability-Probability Plot}
\alias{pp.plot}

\description{
Plots the empirical distribution of a data set against the best fitting
von Mises distribution function.
}

\usage{
pp.plot(x, ref.line = TRUE, tol=1e-20,  xlab = "von Mises Distribution", 
  ylab = "Empirical Distribution", control.circular = list(), ...)
}

\arguments{
\item{x}{a vector. The object is coerced to class \code{\link{circular}}.}
\item{ref.line}{logical, if TRUE a 45 degree reference line is added to
  the plot. Default is TRUE.}
\item{tol}{parameter passed to \code{\link{pvonmises}}.}
\item{xlab, ylab}{labels of the axis.}
\item{control.circular}{the attribute of the resulting object.}
\item{...}{parameters passed to the \code{\link{plot.default}} function.}
}

\value{
a list with the estimated mean and concentration parameter for a von
Mises distribution.
}


\details{
The maximum likelihood estimates of the parameters of the von Mises distribution are computed from the given data set.  The empirical distribution function is plotted against a von Mises distribution function with parameters given by the MLEs computed. 
}

\author{Claudio Agostinelli and Ulric Lund}

\references{
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 10.2, World Scientific Press, Singapore.
}

\seealso{
\code{\link{mle.vonmises}}
}

\examples{
x <- rvonmises(n=25, mu=circular(0), kappa=3)
pp.plot(x)
x <- c(rvonmises(n=20, mu=circular(0), kappa=7), 
  rvonmises(n=20, mu=circular(pi), kappa=7))
pp.plot(x)
}

\keyword{hplot}
