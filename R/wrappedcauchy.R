
###############################################################
#                                                             #
#       Original Splus: Ulric Lund                            #
#       E-mail: ulund@calpoly.edu                             #
#                                                             #
###############################################################

#############################################################
#                                                           #
#   rwrappedcauchy function                                 #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: July, 24, 2003                                    #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1                                             #
#############################################################

rwrappedcauchy <- function(n, mu = 0, rho = exp(-1), units=c("radians", "degrees"), ...) {
    units <- match.arg(units)
    if (units=="degrees") {
        mu <- mu/180*pi
    }
    
    if (rho == 0)
    result <- runif(n, 0, 2 * pi)
    else if (rho == 1)
         result <- rep(mu, n)
    else {
       scale <-  - log(rho)
       result <- rcauchy(n, mu, scale) %% (2 * pi)
    }
    if (units=="degrees") result <- result/pi*180
    result <- circular(result, units=units, ...)
    return(result)
}

#############################################################
#                                                           #
#   dwrappedcauchy function                                 #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: July, 24, 2003                                    #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1                                             #
#############################################################

dwrappedcauchy <- function(x, mu=0, rho=exp(-1)) {
    if (rho < 0 | rho > 1)
        stop("rho must be between 0 and 1")
   
    x <- as.circular(x)
    xcircularp <- circularp(x)
    units <- xcircularp$units
    x <- conversion.circular(x, units="radians")
    n <- length(x)
    attr(x, "circularp") <-  NULL
    
    if (units=="degrees") {
        mu <- mu/180*pi
    }  
    d <- (1 - rho^2)/((2 * pi) * (1 + rho^2 - 2 * rho * cos(x - mu)))
    d <- unclass(d)
    return(d)

  }
