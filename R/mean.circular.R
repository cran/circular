
#############################################################
#                                                           #
#       Original Splus: Ulric Lund                          #
#       E-mail: ulund@calpoly.edu                           #
#                                                           #
#############################################################

#############################################################
#                                                           #
#   mean.circular function                                  #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: August, 10, 2006                                  #
#   Version: 0.4-3                                          #
#                                                           #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#############################################################

mean.circular <- function(x, na.rm=FALSE, control.circular=list(), ...) {
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
   circmean <- MeanCircularRad(x)
   circmean <- conversion.circular(circular(circmean), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
   return(circmean)
}

MeanCircularRad <- function(x) {
   if (any(is.na(x))) {
       circmean <- NA
   } else {
       sinr <- sum(sin(x))
       cosr <- sum(cos(x))
       if (sqrt((sinr^2 + cosr^2))/length(x) > .Machine$double.eps) {
           circmean <- atan2(sinr, cosr)
       } else {
           circmean <- NA
       }
   }
   return(circmean)
}
