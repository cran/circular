### This is necessary since stats::median do not have ... argument
### Work around suggested by Kurt but does not work.
##median <- function(x, na.rm, ...) UseMethod("median")
##median.default <- function(x, na.rm, ...) stats::median(x, na.rm)

#############################################################
#                                                           
#   medianCircular function                                  
#   Author: Claudio Agostinelli                             
#   E-mail: claudio@unive.it                                
#   Date: October, 12, 2009                                  
#   Version: 0.2                                          
#                                                           
#   Copyright (C) 2009 Claudio Agostinelli                  
#                                                           
#############################################################

medianCircular <- function(x, na.rm=FALSE, type="Fisher", deviation=FALSE, control.circular=list(), ...) {
  ## For now only the definition in
  ## equations 2.32 & 2.33
  ## from N.I. Fisher's 'Statistical Analysis of Circular Data',
  ## Cambridge Univ. Press 1993.
  ## is implemented
   type <- match.arg(type)  
   if (na.rm)
       x <- x[!is.na(x)]
   if (length(x)==0) {
        warning("No observations (at least after removing missing values)")
        return(NULL)
   }

   if (is.circular(x)) {
      datacircularp <- circularp(x)
   } else {
      datacircularp <- list(type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="counter")
   }
   dc <- control.circular
   if (is.null(dc$type))
      dc$type <- datacircularp$type
   if (is.null(dc$units))
      dc$units <- datacircularp$units
   if (is.null(dc$template))
      dc$template <- datacircularp$template
   if (is.null(dc$modulo))
      dc$modulo <- datacircularp$modulo
   if (is.null(dc$zero))
      dc$zero <- datacircularp$zero
   if (is.null(dc$rotation))
      dc$rotation <- datacircularp$rotation
   x <- conversion.circular(x, units="radians", zero=0, rotation="counter")
   attr(x, "class") <- attr(x, "circularp") <-  NULL
   if (type=="Fisher")
     circmedian <- MedianFisherCircularRad(x)
   else
     stop("Others 'type' not yet implemented")
   circmedian$median <- conversion.circular(circular(circmedian$median), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
   if (deviation)
     return(circmedian)
   else
     return(circmedian$median)
}

MedianFisherCircularRad <- function(x) {
  dev <- function(x, theta, shift=0) {
  ## x = median
    n <- length(theta)
    res <- pi - sum(abs(pi-abs(MinusPiPlusPiRad(theta-x-shift))))/n
    return(res)
  }
  grid <- res <- seq(0, 2*pi, pi/50)
  for (i in 1:length(grid)) {
    res[i] <- dev(x=grid[i], theta=x)
  }
  
  pos <- which.min(res)
  res <- optim(par=0, fn=dev, lower=MinusPiPlusPiRad(grid[(pos-4)%%length(grid)]-grid[pos]), upper=MinusPiPlusPiRad(grid[(pos+4)%%length(grid)]-grid[pos]), theta=x, shift=grid[pos], method="L-BFGS-B")
  median <- (res$par+grid[pos])%%(2*pi)
  md <- dev(x=median, theta=x, shift=0)
  res <- list(median=median, deviation=md)
  return(res)
}


