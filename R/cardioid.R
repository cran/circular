
###############################################################
#                                                             #
#       Original Splus: Ulric Lund                            #
#       E-mail: ulund@calpoly.edu                             #
#                                                             #
###############################################################

#############################################################
#                                                           #
#   rcardioid function                                      #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: April, 29, 2003                                   #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1                                             #
#############################################################

rcardioid <- function(n, mu=0, rho=0, units=c("radians", "degrees"), ...) {
    units <- match.arg(units)
    if (units=="degrees") {
        mu <- mu/180*pi
    }

    if (rho < -0.5 | rho > 0.5)
        stop("rho must be between -0.5 and 0.5")        
    i <- 1
    result <- rep(0, n)
    while(i <= n) {
        x <- runif(1, 0, 2 * pi)
        y <- runif(1, 0, (1 + 2 * rho)/(2 * pi))
        f <- (1 + 2 * rho * cos(x - mu))/(2 * pi)
        if(y <= f) {
            result[i] <- x
            i <- i + 1
        }
    }    
    if (units=="degrees") result <- result/pi*180
    result <- circular(result, units=units, ...)
    return(result)
}

#############################################################
#                                                           #
#   dcardioid function                                      #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: July, 23, 2003                                    #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1-1                                           #
#############################################################

dcardioid <- function(x, mu=0, rho=0) {
  
    x <- as.circular(x)
    xcircularp <- circularp(x)
    units <- xcircularp$units
    x <- conversion.circular(x, units="radians")
    n <- length(x)
    attr(x, "circularp") <-  NULL
    
    if (units=="degrees") {
        mu <- mu/180*pi
    }  

    if (rho < -0.5 | rho > 0.5)
        stop("rho must be between -0.5 and 0.5")
    d <- (1 + 2 * rho * cos(x - mu))/(2 * pi)
    d <- unclass(d)
    return(d)
}
