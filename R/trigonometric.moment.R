
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
#   Date: July, 24, 2003                                    #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1                                             #
#############################################################

trigonometric.moment <- function(x, p = 1, center = FALSE) {
    x <- as.circular(x)
    xcircularp <- circularp(x)
    units <- xcircularp$units
    x <- conversion.circular(x, units="radians")
    n <- length(x)
    attr(x, "circularp") <-  NULL

    sinr <- sum(sin(x))
    cosr <- sum(cos(x))
    circmean <- atan(sinr, cosr)
    sin.p <- sum(sin(p * (x - circmean * center)))/n
    cos.p <- sum(cos(p * (x - circmean * center)))/n
    mu.p <- atan(sin.p, cos.p)
    rho.p <- sqrt(sin.p^2 + cos.p^2)
    if (units=="degrees") mu.p <- mu.p/pi*180
    attr(mu.p, "circularp") <- xcircularp
    attr(mu.p, "class") <- "circular"
    return(list(mu=mu.p, rho=rho.p, cos=cos.p, sin=sin.p, p=p, n=n, call=match.call()))
}
