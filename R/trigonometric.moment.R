
###############################################################
#                                                             #
#       Original Splus: Ulric Lund                            #
#       E-mail: ulund@calpoly.edu                             #
#                                                             #
###############################################################

#############################################################
#                                                           #
#   trigonometric.moment function                           #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: December, 6, 2005                                    #
#   Copyright (C) 2005 Claudio Agostinelli                  #
#                                                           #
#   Version 0.2-1                                             #
#############################################################

trigonometric.moment <- function(x, p = 1, center = FALSE) {
    x <- unlist(x)
    # Handling missing values
    x <- na.omit(x)
    if ((n <- length(x))==0) {
        warning("No observations (at least after removing missing values)")
        return(NULL)
    }      
    x <- as.circular(x)
    xcircularp <- circularp(x)
    units <- xcircularp$units
    x <- conversion.circular(x, units="radians")
    attr(x, "circularp") <-  NULL

    sinr <- sum(sin(x))
    cosr <- sum(cos(x))
    circmean <- atan2(sinr, cosr)
    sin.p <- sum(sin(p * (x - circmean * center)))/n
    cos.p <- sum(cos(p * (x - circmean * center)))/n
    mu.p <- atan2(sin.p, cos.p)
    rho.p <- sqrt(sin.p^2 + cos.p^2)
    if (units=="degrees") mu.p <- mu.p/pi*180
    attr(mu.p, "circularp") <- xcircularp
    attr(mu.p, "class") <- "circular"
    return(list(mu=mu.p, rho=rho.p, cos=cos.p, sin=sin.p, p=p, n=n, call=match.call()))
}
