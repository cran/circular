\name{plot.edf}
\title{Plot Circular Empirical Distribution Function}
\alias{plot.edf}
\alias{lines.edf}

\description{
Plots the empirical distribution function of a circular data set. 
}

\usage{
\method{plot}{edf}(x, type = "s", xlim = c(0, 2 * pi), ylim = c(0, 1), \dots)
\method{lines}{edf}(x, type = "s", \dots)
}

\arguments{
  \item{x}{vector of circular data measured.}
  \item{type, xlim, ylim}{plotting parameters with useful defaults.}
\item{\dots}{optional graphical parameters.  See help section on \code{\link{par}}.}
}

\note{
Creates a plot or adds a plot (\code{lines.edf}) of the empirical
distribution function of the circular data vector.
}

\details{
The vector of data is taken modulo 2*pi, and then the linear ranks are used to generate an empirical distribution function.
}

\author{Claudio Agostinelli and Ulric Lund}

\seealso{
\code{\link{par}}
}

\examples{
# Compare the edf's of two simulated sets of data.
data1 <- rvonmises(n=10, mu=0, kappa=3)
data2 <- rvonmises(n=10, mu=0, kappa=1)
plot.edf(data1, xlab="Data", ylab="EDF", main="Plots of Two EDF's")
lines.edf(data2, lty=2, col=2)
}

\keyword{hplot}
