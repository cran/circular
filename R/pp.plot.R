
###############################################################
#                                                             #
#       Original Splus: Ulric Lund                            #
#       E-mail: ulund@calpoly.edu                             #
#                                                             #
###############################################################

#############################################################
#                                                           #
#   pp.plot function                                        #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: July, 24, 2003                                    #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1                                             #
#############################################################

pp.plot <- function(x, ref.line = TRUE, tol=1e-20,  xlab = "von Mises Distribution", ylab = "Empirical Distribution", ...) {
        x <- as.circular(x)
        xcircularp <- circularp(x)
        units <- xcircularp$units
        x <- conversion.circular(x, units="radians")

        res <- mle.vonmises(x)
        mu <- res$mu
    kappa <- res$kappa

    n <- length(x)
    x <- sort(x %% (2 * pi))
    z <- (1:n)/(n + 1)
    
    y <- pvonmises(q=x, mu=mu, kappa=kappa, tol=tol)
    
    plot.default(z, y, xlab=xlab, ylab=ylab, ...)
    if (ref.line)
        abline(0, 1)
        if (units=="degrees") mu <- mu/pi*180
        attr(mu, "circularp") <- xcircularp
        attr(mu, "class") <- "circular"
    invisible(list(mu=mu, kappa=kappa))
}
