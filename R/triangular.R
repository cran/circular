
###############################################################
#                                                             #
#       Original Splus: Ulric Lund                            #
#       E-mail: ulund@calpoly.edu                             #
#                                                             #
###############################################################

#############################################################
#                                                           #
#   rtriangular function                                    #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: July, 24, 2003                                    #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1                                             #
#############################################################

rtriangular <- function(n, rho, units=c("radians", "degrees"), ...) {
        if (rho < 0 | rho > 4/pi^2) stop("'rho' must be between 0 and 4/pi^2")
        units <- match.arg(units)
    u <- matrix(c(runif(n)), ncol = 1)
    get.theta <- function(u, rho)
    {
        if(u < 0.5) {
            a <- pi * rho
            b <-  - (4 + pi^2 * rho)
            c <- 8 * pi * u
            theta1 <- ( - b + sqrt(b^2 - 4 * a * c))/(2 * a)
            theta2 <- ( - b - sqrt(b^2 - 4 * a * c))/(2 * a)
            min(theta1, theta2)
        }
        else {
            a <- pi * rho
            b <- 4 - 3 * pi^2 * rho
            c <- (2 * pi^3 * rho) - (8 * pi * u)
            theta1 <- ( - b + sqrt(b^2 - 4 * a * c))/(2 * a)
            theta2 <- ( - b - sqrt(b^2 - 4 * a * c))/(2 * a)
            max(theta1, theta2)
        }
    }
    theta <- apply(u, 1, get.theta, rho)
    theta[theta > pi] <- theta[theta > pi] - 2 * pi
        if (units=="degrees") theta <- theta*180/pi
    return(circular(theta, units=units, ...))
}

#############################################################
#                                                           #
#   dtriangular function                                    #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: July, 24, 2003                                    #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1                                             #
#############################################################

dtriangular <- function(x, rho) {
    if (rho < 0 | rho > 4/pi^2) stop("'rho' must be between 0 and 4/pi^2")
    x <- as.circular(x)
    xcircularp <- circularp(x)
    units <- xcircularp$units
    x <- conversion.circular(x, units="radians")
    n <- length(x)
    attr(x, "circularp") <-  NULL
    d <- (4 - pi^2 * rho + 2 * pi * rho * abs(pi - x))/(8 * pi)
    d <- unclass(d)
    return(d)

}
