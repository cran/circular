
###############################################################
#                                                             #
#       Original Splus: Ulric Lund                            #
#       E-mail: ulund@calpoly.edu                             #
#                                                             #
###############################################################

#############################################################
#                                                           #
#   mle.vonmises function                                   #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: September, 22, 2003                               #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#   Version 0.2-3                                           #
#############################################################

mle.vonmises <- function(x, mu, kappa, bias=FALSE) {

    x <- as.circular(x)
    xcircularp <- circularp(x)
    units <- xcircularp$units
    x <- conversion.circular(x, units="radians")

    n <- length(x)
    sinr <- sum(sin(x))
    cosr <- sum(cos(x))
    est.mu <- FALSE 
    if (missing(mu)) {  
        mu <- atan(sinr, cosr)
        est.mu <- TRUE
    } else {
        if (units=="degrees") mu <- mu/180*pi
    }
    est.kappa <- FALSE
    if (missing(kappa)) {
        V <- mean.default(cos(x - mu))
        if (V > 0) {
            kappa <- A1inv(V)
        } else {
            kappa <- 0
        }
        if (bias == TRUE) {
            if (kappa < 2) {
                kappa <- max(kappa - 2 * (n * kappa)^-1, 0)
            } else {
                kappa <- ((n - 1)^3 * kappa)/(n^3 + n)
            }
        }
        est.kappa <- TRUE
    }

    A1temp <- A1(kappa)
    se.mu <- se.kappa <- 0
    if (est.mu) se.mu <- 1/(n*kappa*A1temp)
    if (est.kappa) se.kappa <- 1/(n*(1-A1temp/kappa-A1temp^2))
    result <- list()

    if (units=="degrees") {
        mu <- mu/pi*180
    }
    
    attr(mu, "circularp") <- xcircularp
    attr(mu, "class") <- "circular"
    
    result$call <- match.call()
    result$mu <- mu
    result$kappa <- kappa
    result$se.mu <- sqrt(se.mu)
    result$se.kappa <- sqrt(se.kappa)
    result$est.mu <- est.mu
    result$est.kappa <- est.kappa
    class(result) <- "mle.vonmises"
    return(result)
}

#############################################################
#                                                           #
#   print.mle.vonmises function                             #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: November, 19, 2003                                #
#   Version: 0.1-2                                          #
#                                                           #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################

print.mle.vonmises <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    cat("mu: ")
    cat(format(x$mu, digits=digits), " (", format(x$se.mu, digits=digits), ")\n")
    cat("\n")
    cat("kappa: ")    
    cat(format(x$kappa, digits=digits), " (", format(x$se.kappa, digits=digits), ")\n")
    cat("\n")    
    if (!x$est.mu) cat("mu is known\n")
    if (!x$est.kappa) cat("kappa is known\n")
    invisible(x)
}
