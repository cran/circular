#############################################################
#                                                           #
#   as.circular function                                    #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: December, 5, 2005                                 #
#   Copyright (C) 2005 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1-1                                           #
#############################################################

as.circular <- function (x, ...) {
    if (is.circular(x)) return(x)
    else if(!is.null(xcircularp <- circularp(x))) circular(x, type=xcircularp$type, units=xcircularp$units, template=xcircularp$template, modulo=xcircularp$modulo, zero=xcircularp$zero, rotation=xcircularp$rotation)
    else {
           if (is.null(list(..., expand.dots=TRUE)$units)) 
               warning("object 'x' is coerced to the class 'circular' using radian unit of measure in: ", sys.call(-1))
           circular(x, ...)
         }
}

