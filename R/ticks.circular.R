#############################################################
#                                                           #
#   ticks.circular function                                 #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: May, 29, 2006                                     #
#   Version: 0.3-1                                            #
#                                                           #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#############################################################
 
ticks.circular <- function(x, template=c("none", "geographics"), zero=NULL, rotation=NULL, tcl=0.025, col=NULL, ...) {
  template <- match.arg(template)
  if (is.null(col)) col <- par("col")
  x <- conversion.circular(x, units="radians", template=template, zero=zero, rotation=rotation, modulo="2pi")
  attr(x, "circularp") <- attr(x, "class") <- NULL
  TicksCircularRad(x, tcl, col, ...)
}

TicksCircularRad <- function(x, tcl, col, ...) {
    r <- 1+tcl*c(-1/2,1/2)
    z <- cos(x)
    y <- sin(x)
    for (i in 1:length(x)) {
       lines.default(z[i]*r, y[i]*r, col=col, ...)
    }
}

