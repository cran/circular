
###############################################################
#                                                             #
#       Original Splus: Ulric Lund                            #
#       E-mail: ulund@calpoly.edu                             #
#                                                             #
###############################################################

#############################################################
#                                                           #
#   rwrappednormal function                                 #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: December, 17, 2003                                #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1-1                                           #
#############################################################

rwrappednormal <- function(n, mu=0, rho, sd=1, units=c("radians", "degrees"), ...) {
    units <- match.arg(units)
    if (units=="degrees") {
        mu <- mu/180*pi
    }
    if (missing(rho)) {
        rho <- exp(-sd^2/2)
    }
    if (rho < 0 | rho > 1)
        stop("rho must be between 0 and 1")        
    if (rho == 0)
        result <- runif(n, 0, 2 * pi)
    else if (rho == 1)
        result <- rep(mu, n)
    else {
        sd <- sqrt(-2 * log(rho))
        result <- rnorm(n, mu, sd) %% (2 * pi)
    }
    if (units=="degrees") result <- result/pi*180
    result <- circular(result, units=units, ...)
    return(result)
}

#############################################################
#                                                           #
#   dwrappednormal function                                 #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: December, 17, 2003                                #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1-2                                           #
#############################################################

dwrappednormal <- function(x, mu=0, rho, sd=1, K, min.k=10) {
  
    x <- as.circular(x)
    xcircularp <- circularp(x)
    units <- xcircularp$units
    x <- conversion.circular(x, units="radians")
    n <- length(x)
    im <- length(mu)
    attr(x, "class") <- attr(x, "circularp") <-  NULL
    
    if (units=="degrees") {
        mu <- mu/180*pi
    }  
    if (missing(rho)) {
        rho <- exp(-sd^2/2)
    }
    if (rho < 0 | rho > 1)
        stop("rho must be between 0 and 1")
    var <- -2 * log(rho)
    sd <- sqrt(var)
    
    if (missing(K)) {
        range <- abs(mu-x)
        K <- (range+6*sqrt(var))%/%(2*pi)+1
        K <- max(min.k, K)
    }

    z <- .Fortran("dwrpnorm",
                    as.double(x),
                    as.double(mu),
                    as.double(sd),
                    as.integer(n),
                    as.integer(im),
                    as.integer(K),
                    d=mat.or.vec(im, n),
                    PACKAGE="circular"
    )
    d <- z$d/sqrt(var * 2 * pi)

## we have to be careful to return the right dimension when x is a matrix and mu is a vector 
#    if (!is.null(dim(x))) {
#         d <- matrix(d, dim(x))
#    }
    return(d)
}
