
###############################################################
#                                                             #
#       Original Splus: Ulric Lund                            #
#       E-mail: ulund@calpoly.edu                             #
#                                                             #
###############################################################

#############################################################
#                                                           #
#   rho.circular function                                   #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: April, 11, 2005                                   #
#   Version: 0.2                                            #
#                                                           #
#   Copyright (C) 2005 Claudio Agostinelli                  #
#                                                           #
#############################################################

rho.circular <- function(x, na.rm=FALSE) {
    if (na.rm) 
        x <- x[!is.na(x)]
    if (any(is.na(x))) {
        warning("No observations (at least after removing missing values)")
        return(NA)
    }

    n <- length(x)
    x <- as.circular(x)
    x <- conversion.circular(x, units="radians")
    
    sinr <- sum(sin(x))
    cosr <- sum(cos(x))
    result <- sqrt(sinr^2 + cosr^2)/n
    return(result)
}
