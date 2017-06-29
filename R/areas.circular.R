#############################################################
#
#	areas.region.circular
#       GNU General Public Licence 2.0 
#	Author: Claudio Agostinelli
#	E-mail: claudio@unive.it
#	Date: August, 26, 2013
#	Version: 0.1
#
#	Copyright (C) 2013 Claudio Agostinelli
#
#############################################################

areas.region.circular  <- function(x, breaks=NULL, z=NULL, bw, adjust = 1, type = c("K", "L"), kernel = c("vonmises", "wrappednormal"), na.rm = FALSE, step=0.01, ...) {
  if (is.null(z))
    z <- circular(seq(-step,2*pi+step,step))
  if (is.null(breaks))
    breaks <- circular(seq(0, 2*pi+pi/4, pi/4))
  if (is.circular(x))
    xcp <- circularp(x)
  else
    xcp <- list(type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="counter")   
  x <- conversion.circular(x, units="radians", zero=0, rotation="counter", modulo="2pi")
  z <- conversion.circular(z, units="radians", zero=0, rotation="counter", modulo="asis")
  if (is.vector(breaks))
    modulobreaks <- "2pi"
  else
    modulobreaks <- "asis"
  breaks <- conversion.circular(breaks, units="radians", zero=0, rotation="counter", modulo=modulobreaks)    
  class(breaks) <- class(breaks)[class(breaks)!="circular"]
  attr(breaks, "circularp") <- NULL
  if (is.vector(breaks)) {
    breaks <- sort(unique(breaks))
    breaks <- cbind(breaks, c(breaks[-1], breaks[1]+2*pi))
  }
  extend <- range(breaks)
  extend <- c(floor(extend[1]/(2*pi)), ceiling(extend[2]/(2*pi)))
  object <- density.circular(x=x, z=z, bw=bw, adjust=adjust, type=type, kernel=kernel, na.rm=na.rm)
  areas <- area2(breaks, object, extend=extend)
  result <- list()
  xunits <- circularp(x)$units
  result$breaks <- conversion.circular(circular(breaks), xcp$units, xcp$type, xcp$template, 'asis', xcp$zero, xcp$rotation)
  result$tot <- areas$tot
  result$areas <- areas$areas
  object$x <- conversion.circular(object$x, xcp$units, xcp$type, xcp$template, 'asis', xcp$zero, xcp$rotation)
  result$density <- object
  class(result) <- 'areas.region.circular'
  return(result)
}

## Calculate areas under several disjoint intervals
area2 <- function(x, object, extend, ...) {
#x: is a matrix with two columns
#object: an object from density.circular
#...: values passed to integrate function
  extend <- seq(extend[1], extend[2],1)
###  extend <- extend[extend!=0]
  den <- approxfun(x=as.vector(outer(object$x,extend*2*pi,FUN="+")), y=rep(object$y, length(extend)))
  int <- function(x) integrate(f=den, lower=x[1], upper=x[2], ...)$value
  areas <- apply(x, 1, int)
  tot <- sum(areas)
  result <- list(tot=tot, areas=areas)
  return(result)
}

if (FALSE) {
#### EXAMPLES
set.seed(1234)
x <- c(rvonmises(100, circular(0), 8, control.circular=list(units="hours")), rvonmises(100, circular(pi), 8, control.circular=list(units="hours")))
plot(x, template="clock24")
res1 <- areas.region.circular(x, breaks=circular(c(5,8,17,19,21), units="hours"), bw=10) ## a partition of the circle
res2 <- areas.region.circular(x, breaks=circular(matrix(c(7,18), nrow=1), units="hours"), bw=10) ## a single interval 
res3 <- areas.region.circular(x, breaks=circular(matrix(c(7,18,19,6+24), ncol=2), units="hours"), bw=10) ## a two or more intervals
res4 <- areas.region.circular(x, breaks=circular(matrix(c(6,6+3*24), ncol=2), units="hours"), bw=10) ## over more than one clock
res5 <- areas.region.circular(x, breaks=circular(matrix(c(6,6+3*24), ncol=2), units="hours"), bw=10, step=0.0001) ## increase precision as in modal.region.circular
}
