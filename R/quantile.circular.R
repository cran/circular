#############################################################
#                                                           #
#   quantile.circular function                              #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: July, 18, 2011                                    #
#   Copyright (C) 2011 Claudio Agostinelli                  #
#                                                           #
#   Version 0.2                                             #
#############################################################

quantile.circular <- function(x, probs = seq(0, 1, 0.25), na.rm = FALSE, names = TRUE, type = 7, control.circular=list(), ...) {
# Handling missing values
  x <- na.omit(x)
  if (length(x)==0) {
    warning("No observations (at least after removing missing values)")
    return(NULL)
  }
  if (is.circular(x))
    datacircularp <- circularp(x)     
  else
    datacircularp <- list(type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="counter")   
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
  x <- conversion.circular(x, units="radians", zero=0, rotation="counter", modulo="2pi")
  attr(x, "class") <- attr(x, "circularp") <- NULL
  q <- QuantileCircularRad(x=x, probs=probs, names=names, type=type)
  q <- conversion.circular(circular(q), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
  return(q)
}

QuantileCircularRad <- function(x, probs, names, type) {
  n <- length(x)
  x <- sort(x %% (2 * pi))
  spacings <- c(diff(x), x[1L] - x[n] + 2 * pi)
  max.spacing <- (1:n)[spacings == max(spacings)]
  off.set <- 2 * pi - x[max.spacing + 1]
  if (max.spacing != n)
    x2 <- x + off.set
  else
    x2 <- x
  x2 <- sort(x2 %% (2 * pi))
  q <- quantile.default(x=x2, probs=probs, names=names, type=type)
  if (max.spacing != n)
    q <- q - off.set
  return(q)
}
