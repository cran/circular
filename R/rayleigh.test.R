
###############################################################
#                                                             #
#       Original Splus: Ulric Lund                            #
#       E-mail: ulund@calpoly.edu                             #
#                                                             #
###############################################################

#############################################################
#                                                           #
#   rayleigh.test function                                  #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: July, 25, 2003                                    #
#   Version: 0.1                                            #
#                                                           #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################

rayleigh.test <- function(x, mu) {
    x <- as.circular(x)
    xcircularp <- circularp(x)
    units <- xcircularp$units
    x <- conversion.circular(x, units="radians")
    attr(x, "circularp") <- attr(x, "class") <- NULL
    n <- length(x)
    if (missing(mu)) {
        ss <- sum(sin(x))
        cc <- sum(cos(x))
        rbar <- (sqrt(ss^2 + cc^2))/n
        z <- (n * rbar^2)
        p.value <- exp( - z)
        if (n < 50)
        temp <- (1 + (2 * z - z^2)/(4 * n) - (24 * z - 132 * z^2 + 76 * z^3 - 9 * z^4)/(288 * n^2))
        else temp <- 1
        result <- list(statistic = rbar, p.value = p.value * temp, mu=NA)
    } else {
    r0.bar <- (sum(cos(x - mu)))/n
    z0 <- sqrt(2 * n) * r0.bar
    pz <- pnorm(z0)
    fz <- dnorm(z0)
    p.value <- 1 - pz + fz * ((3 * z0 - z0^3)/(16 * n) + (15 * z0 + 305 * z0^3 - 125 * z0^5 + 9 * z0^7)/(4608 * n^2))
    result <- list(statistic = r0.bar, p.value = p.value, mu=mu)
    }
    class(result) <- "rayleigh.test"
    return(result)
}

#############################################################
#                                                           #
#   print.rayleigh.test function                            #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: November, 19, 2003                                #
#   Version: 0.1-1                                          #
#                                                           #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################

print.rayleigh.test <- function(x, digits=4, ...) {
    rbar <- x$statistic
    p.value <- x$p.value
    mu <- x$mu
    cat("\n", "      Rayleigh Test of Uniformity \n")
    if (is.na(mu)) {
        cat("       General Unimodal Alternative \n\n")
    } else {
        cat("       Alternative with Specified Mean Direction: ", mu, "\n\n")
    }
    cat("Test Statistic: ", round(rbar, digits=digits), "\n")
    cat("P-value: ", round(p.value, digits=digits), "\n\n")
    invisible(x)
}
