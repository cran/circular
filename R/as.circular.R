#############################################################
#                                                           #
#   as.circular function                                    #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: July, 24, 2003                                    #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1                                             #
#############################################################

as.circular <- function (x, ...) {
    if (is.circular(x)) return(x)
    else if(!is.null(xcircularp <- circularp(x))) circular(x, type=xcircularp$type, units=xcircularp$units, template=xcircularp$template, modulo=xcircularp$modulo, zero=xcircularp$zero, rotation=xcircularp$rotation)
    else circular(x, ...)
}

