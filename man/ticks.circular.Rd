\name{ticks.circular}
\alias{ticks.circular}
\title{Draw Tick-Marks in a Circular Plot}
\description{
Draw tick-marks in a circular plot.
}
\usage{
ticks.circular(x, template = c("none", "geographics"), zero, rotation, tcl = 0.025, col, \dots)
}

\arguments{

  \item{x}{the points at which tick-marks are to be drawn.}
\item{template}{either \code{none} or \code{geographics}.}
\item{zero}{the zero of the plot (in radians).}
\item{rotation}{the rotation of the plot.}
\item{col}{color for the tick marks. If missing, function uses `par("col.axis")'. }
\item{tcl}{The length of tick marks.}
\item{\dots}{parameters passed to \code{line.default}.}
}

\author{Claudio Agostinelli}

\seealso{
  \code{\link{plot.circular}} and \code{\link{axis.circular}}
}

\keyword{hplot}
