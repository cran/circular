
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
#   Date: July, 24, 2003                                    #
#   Version: 0.1                                            #
#                                                           #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################

rho.circular <- function(x, na.rm=FALSE) {
    if (na.rm) 
        x <- x[!is.na(x)]
    if (any(is.na(x))) return(NA)

    n <- length(x)
    x <- as.circular(x)
    x <- conversion.circular(x, units="radians")
    
    sinr <- sum(sin(x))
    cosr <- sum(cos(x))
    sqrt(sinr^2 + cosr^2)/n
}
