#############################################################
#  weighted.mean.circular function
#  Author: Claudio Agostinelli
#  E-mail: claudio@unive.it
#  Date: May, 12, 2015
#  Version: 0.1
#  Copyright (C) 2015 Claudio Agostinelli
#############################################################

weighted.mean.circular <- function(x, w, na.rm=FALSE, control.circular=list(), ...) {
  if (missing(w))
    mean.circular(x=x, na.rm=na.rm, control.circular=control.circular, ...)
  if (any(is.na(w))) {
    warning("Missing values are not allowed in the weights vector")
    return(circular(NA))
  }
  if (length(w) != length(x)) 
    stop("'x' and 'w' must have the same length")
  w <- as.double(w)
  if (na.rm) {
    nax <- !is.na(x)
    x <- x[nax]
    w <- w[nax]
  }
  neq0 <- w !=0
  x <- x[neq0]
  w <- w[neq0]
   
  if (length(x)==0) {
    warning("No observations (at least after removing missing values)")
    return(circular(NA))
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
  x <- conversion.circular(x, units="radians")
  attr(x, "class") <- attr(x, "circularp") <-  NULL
  circmean <- WeightedMeanCircularRad(x, w)
  circmean <- conversion.circular(circular(circmean, template=datacircularp$template, zero=datacircularp$zero, rotation=datacircularp$rotation), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
  return(circmean)
}

WeightedMeanCircularRad <- function(x, w) {
  if (any(is.na(x))) {
    circmean <- NA
  } else {
    circmean <- .C("WeightedMeanCircularRad",
      x=as.double(x),
      w=as.double(w),
      n=as.integer(length(x)),
      result=as.double(0),
      PACKAGE="circular")$result
  }
  return(circmean)
}
