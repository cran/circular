#############################################################
#                                                           #
#       Original Splus: Ulric Lund                          #
#       E-mail: ulund@calpoly.edu                           #
#                                                           #
#############################################################


#############################################################
#                                                           #
#   range.circular function                                 #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: August, 10, 2006                                  #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#   Version 0.3-2                                           #
#############################################################

range.circular <- function(x, test = FALSE, na.rm=FALSE, finite=FALSE, control.circular=list(), ...) {
   if (finite) 
       x <- x[is.finite(x)]
   else if (na.rm) 
       x <- x[!is.na(x)]

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
   
    x <- conversion.circular(x, units="radians", zero=0, rotation="counter", modulo="2pi")
    attr(x, "class") <- attr(x, "circularp") <- NULL

    result <- RangeCircularRad(x, test)
   
    if (test) {
       result$range <- conversion.circular(circular(result$range), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
    } else {
      result <- conversion.circular(circular(result), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
    }
    return(result)
}

RangeCircularRad <- function(x, test=TRUE) {
   x <- sort(x %% (2*pi))
   n <- length(x)
   spacings <- c(diff(x), x[1] - x[n] + 2*pi)
   range <- 2*pi - max(spacings)
   if (test == TRUE) {
       stop <- floor(1/(1 - range/(2*pi)))
       index <- c(1:stop)
       sequence <- ((-1)^(index - 1)) * exp(log(gamma(n + 1)) - log(gamma(index + 1)) - log(gamma(n - index + 1))) * (1 - index * (1 - range/(2 * pi)))^(n - 1)
       p.value <- sum(sequence)
       result <- list(range=range, p.value=p.value)
   } else {
       result <- range
   }
   return(result)
}
