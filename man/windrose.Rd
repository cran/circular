\name{windrose}
\alias{windrose}
\title{Windrose Generator}
       \description{This function creates a windrose used to visualize
	 the direction and magnitude of wind.  The pedals of a windrose
	 indicate the proportion of time wind comes from a given
	 direction. Bands on the windrose indicate the proportions of
	 winds of each magnitude.}
\usage{windrose(x, y=NULL, breaks=NULL, bins=12, increment = 10, 
  main='Wind Rose', cir.ind = 0.05, fill.col=NULL, plot.mids=TRUE, 
  mids.size=1.2, osize=0.1, axes=TRUE, ticks=TRUE, tcl=0.025, 
  tcl.text=-0.15, cex = 1, digits=2, num.ticks=12, xlim=c(-1.2, 1.2), 
  ylim=c(-1.2, 1.2), uin, tol=0.04, right=FALSE, shrink=NULL, 
  label.freq=FALSE, calm=c("0", "NA"), ...)
}

\arguments{
         \item{x}{a vector contains direction or a two columns data frame, where the first component is the direction and the second the magnitude.  The vector or the first column in the case of data frame is coerced to class \code{\link{circular}}.}
         \item{y}{a vector contains magnitude. If 'y' is not NULL and 'x' is a dataframe,  only the first column of 'x' is used for direction.}
         \item{breaks}{the extremes of the pedals. The biggest value (in 2*pi) is recycled for building the first pedal. The vector is coerced to class \code{\link{circular}} but only the units is used.}
	 \item{bins}{Number of pedals. Ignored if 'breaks' is not NULL.}
	 \item{increment}{Grouping size of magnitude. These are the bins
	 of the magnitudes displayed on each pedal.}
	 \item{main}{Title for plot.}
	 \item{cir.ind}{Percent intervals expressed on each circle if the pedals are equally spaced, otherwise values of density}
	 \item{fill.col}{colors used to fill the pedals for each magnitude. The colors are recycled if necessary. The default is to use 'blue' and 'red'.} 
         \item{plot.mids}{plot lines at the midpoints of the pedals.}
         \item{mids.size}{length of the lines for midpoints.}
         \item{osize}{radius of the circle draws at the center of the plot.}
         \item{axes}{if TRUE axes are added to the plot. The function \code{\link{axis.circular}} is used.}
         \item{ticks}{if TRUE ticks are added to the plot. The function \code{\link{ticks.circular}} is used.}
         \item{tcl}{length of the ticks.}
         \item{tcl.text}{The position of the axis labels.}
         \item{cex}{point character size.  See help on \code{\link{par}}.}
         \item{digits}{number of digits used to print axis values and other numbers.}
         \item{num.ticks}{number of tick marks draw.}
         \item{tol}{proportion of white space at the margins of plot}
         \item{uin}{desired values for the units per inch parameter. If of length 1, the desired units per inch on the x axis.}
         \item{xlim, ylim}{the ranges to be encompassed by the x and y axes. Useful for centering the plot.}
         \item{right}{logical; if TRUE, the pedals are right-closed (left open) intervals.}
         \item{shrink}{maximum length of the pedals, it can be used to plot several graphics with the same scale.}
         \item{label.freq}{logical; if TRUE, the relative frequencies are used in the magnitude instead of intensities, when the breaks are equally spaced.}
         \item{calm}{"0" or "NA", see details below.}
         \item{\dots}{further parameters ignored for now.}
}

\details{Following the convention of the National Weather Service, winds
  with a direction of 0 are considered calm, while winds with a
  direction of 360 degrees (2*pi radians) are assumed to be from the north.  Calm winds
  are excluded from the wind rose creation. We allow, in direction, to use \code{NA} to indicate calm wind (argument \code{calm}). 

  This wind rose preserve areas of pedals, that is counts are proportional to the area of the pedals rather than to the length of the pedals. This is also for the slides created for the magnitudes.
}

\value{
         \item{x}{directions}
         \item{y}{magnitudes}
  	 \item{table}{Matrix output of the counts of wind direction and magnitude.
  	   Columns are in the same units as the data, according to step size, and rows are
  	   based on the increment size.}
	 \item{number.obs}{Total number of observations.}
	 \item{number.calm}{The number of calm observations omitted from the wind rose plot.}
         \item{breaks}{extremes of the pedals.}
         \item{mids}{midpoints of pedals.}
         \item{call}{the \code{\link{match.call}} result.}
}

\note{some codes from \code{\link{eqscplot}} in 'MASS' is used.}
	 
\author{Matt Pocernich <pocernic@rap.ucar.edu>, ported in the package 'circular' by Claudio Agostinelli}

\examples{
# Random distribution of direction and magnitude in degrees

dir <- circular(runif(100, 0, 360), units="degrees")
mag <-  rgamma(100, 15)
sample <- data.frame(dir=dir, mag=mag)

par(mfrow=c(2,2))
res <- windrose(sample)
## we join two pedals and keep the same shrink (scale of the plot)
breaks <-circular(seq(0, 2 * pi, by = pi/6))
breaks <- breaks[-2]
windrose(sample, breaks=breaks, main="The same but with two pedals joined", 
  shrink=res$shrink)
## change the rotation
sample <- data.frame(dir=circular(dir, units="degrees", rotation="clock"), mag=mag)
windrose(sample, breaks=breaks, main="Change the rotation", shrink=res$shrink)
## use geographics template
sample <- data.frame(dir=circular(dir, units="degrees", template="geographics"),
  mag=mag)
windrose(sample, breaks=breaks, main="Use the template 'geographics'", 
  shrink=res$shrink)

## do the same plot but in radians
dir <- conversion.circular(dir)
windrose(x=dir, y=mag, xlim=c(-1.3, 1.3))

## magnify some part of the plot
windrose(x=dir, y=mag, xlim=c(0, 1.3))
}

\keyword{hplot}
