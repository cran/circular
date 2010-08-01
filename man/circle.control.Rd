\name{circle.control}
\alias{circle.control}
\title{Auxiliary for Controlling Circular Plots}
\description{
  Auxiliary function as user interface for circular plots. Typically only used when calling plot.circular. 
}
\usage{
circle.control(n = 1000, type = "l", col = 1, bg = par("bg"), 
  pch = 1, cex = 1, lty = 1, lwd = 1)
}
\arguments{
  \item{n}{number of points used to interpolate the circle}
  \item{type}{1-character string giving the type of plot desired. The following values are possible, for details, see \code{\link{plot}}: "p" for points, "l" for lines, "o" for overplotted points and lines, "b", "c" for (empty if "c") points joined by lines, "s" and "S" for stair steps and "h" for histogram-like vertical lines. Finally, "n" does not produce any points or lines.}
  \item{col}{The color used.}
  \item{bg}{The color to be used for the background of the device region. }
  \item{pch}{Either an integer specifying a symbol or a single character to be used as the default in plotting points. See \code{\link{points}} for possible values and their interpretation. Note that only integers and single-character strings can be set as a graphics parameter (and not \code{NA} nor \code{NULL}).}
  \item{cex}{A numerical value giving the amount by which plotting text and symbols should be magnified relative to the default.}
  \item{lty}{The line type. Line types can either be specified as an integer (0=blank, 1=solid (default), 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash) or as one of the character strings "blank", "solid", "dashed", "dotted", "dotdash", "longdash", or "twodash", where "blank" uses 'invisible lines' (i.e., does not draw them). Alternatively, a string of up to 8 characters (from c(1:9, "A":"F")) may be given, giving the length of line segments which are alternatively drawn and skipped. See section 'Line Type Specification'. }
  \item{lwd}{The line width, a positive number, defaulting to 1. The interpretation is device-specific, and some devices do not implement line widths less than one. (See the help on the device for details of the interpretation.) }
}

\author{Claudio Agostinelli}

\seealso{\code{\link{plot.circular}}}

\examples{
  plot(rvonmises(10, circular(0), 1), control.circle=circle.control(col=2, lty=2))
}

\keyword{hplot}
