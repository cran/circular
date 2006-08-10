\name{axis.circular}
\title{Add Axis to a Circular Plot}
\alias{axis.circular}

\description{
  Add axis to a plot of circular data points on the current graphics
  device.
}

\usage{
axis.circular(at=NULL, labels=NULL,  units = NULL, template=NULL,  
          modulo=NULL, zero=NULL, rotation=NULL, tick=TRUE, lty, lwd, 
          cex, col, font, tcl=0.025, tcl.text=0.125, digits=2)
}

\arguments{
\item{at}{the points at which tick-marks are to be drawn. If \code{NULL} the
tick-marks are placed to 0, pi/2, pi and 3pi/2 radians.}
\item{labels}{a vector of character strings to be placed at the
  tickpoints. If \code{NULL} the labels are choosen according to
  \code{units} and \code{template} arguments.}
\item{units}{either \code{radians} or \code{degrees}. If \code{NULL} the value is taken from the attributes of the object \code{at}.}
\item{template}{either \code{none} or \code{geographics}. If \code{NULL} the value is taken from the attributes of the object \code{at}.}
\item{modulo}{either \code{asis} or \code{2pi} or \code{pi}. If \code{NULL} the value is taken from the attributes of the object \code{at}.}
\item{zero}{the zero of the plot (in radians, counterclockwise). If \code{NULL} the value is taken from the attributes of the object \code{at}.}
\item{rotation}{the rotation of the plot. If \code{NULL} the value is taken from the attributes of the object \code{at}.}
\item{tick}{logical: if \code{TRUE} ticks are plotted at tick-marks.}
\item{lty, lwd}{line type, width for the tick marks. If missing means to
  use `par("lty")' and `par("lwd")'.}
\item{cex}{a numerical value giving the amount by which plotting text and symbols should be scaled relative to the default.}
\item{col}{color for the the tick marks. If missing means to use `par("col.axis")'. }
\item{font}{font for text. If missing means to use `par("font.axis")'.}
\item{tcl}{The length of tick marks.}
\item{tcl.text}{The position of the axis labels.}
\item{digits}{number of digits used to print axis values.}
}

\author{Claudio Agostinelli}

\seealso{
\code{\link{plot.circular}} and \code{\link{ticks.circular}}.
}

\examples{
data.vm <- rvonmises(n=100, mu=circular(0), kappa=3) 
plot(data.vm, axes=FALSE, ticks=FALSE)
axis.circular(at=circular(seq(0, 11/6*pi, pi/6)), labels=c("0",
expression(frac(pi,6)), expression(paste(frac(1,3), pi)),
expression(frac(pi,2)), expression(paste(frac(2,3), pi)),
expression(paste(frac(5,6), pi)), expression(pi),
expression(paste(frac(7,6), pi)), expression(paste(frac(4,3), pi)),
expression(paste(frac(3,2), pi)), expression(paste(frac(5,3), pi)),
expression(paste(frac(11,6), pi))))
}

\keyword{hplot}
