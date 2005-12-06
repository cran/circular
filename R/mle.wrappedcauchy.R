
###############################################################
#                                                             #
#       Original Splus: Ulric Lund                            #
#       E-mail: ulund@calpoly.edu                             #
#                                                             #
###############################################################

#############################################################
#                                                           #
#   mle.wrappedcauchy function                              #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: December, 6, 2005                                   #
#   Copyright (C) 2005 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1-3                                           #
#############################################################

mle.wrappedcauchy <- function(x, mu, rho, tol = 1e-015, max.iter = 100) {

    # Handling missing values
    x <- na.omit(x)
    if (length(x)==0) {
        warning("No observations (at least after removing missing values)")
        return(NULL)
    }

    if (length(tol)==1) tol <- rep(tol, 2)
    if (length(tol) > 2) stop("'tol' must have less than 2 elements")
    x <- as.circular(x)
    xcircularp <- circularp(x)
    units <- xcircularp$units
    x <- conversion.circular(x, units="radians")
    if (missing(mu)) {
        mu <- mean.circular(x)
    } else {
        if (units=="degrees") mu <- mu/180*pi
    }
    
    if (missing(rho)) rho <- rho.circular(x)
    if (rho < 0 | rho > 1) stop("'rho' must be between 0 and 1")
    
    mu1.old <- (2 * rho * cos(mu))/(1 + rho^2)
    mu2.old <- (2 * rho * sin(mu))/(1 + rho^2)
    w.old <- 1/(1 - mu1.old * cos(x) - mu2.old * sin(x))
    flag <- TRUE
    iter <- 0 
    while (flag & iter <= max.iter) {
           iter <- iter + 1
       mu1.new <- sum(w.old * cos(x))/sum(w.old)
       mu2.new <- sum(w.old * sin(x))/sum(w.old)
       diff1 <- abs(mu1.new - mu1.old)
       diff2 <- abs(mu2.new - mu2.old)
       if ((diff1 < tol[1]) && (diff2 < tol[2]))
           flag <- FALSE
       else {
           mu1.old <- mu1.new
           mu2.old <- mu2.new
           w.old <- 1/(1 - mu1.old * cos(x) - mu2.old * sin(x))
       }
    }
    mu.const <- sqrt(mu1.new^2 + mu2.new^2)
    rho <- (1 - sqrt(1 - mu.const^2))/mu.const
    mu <- atan2(mu2.new, mu1.new) %% (2 * pi)
    if (units=="degrees") {
        mu <- mu/pi*180
    }
    
    attr(mu, "circularp") <- xcircularp
    attr(mu, "class") <- "circular"

    result <- list()
    result$call <- match.call()
    result$mu <- mu 
    result$rho <- rho
    result$convergence <- TRUE
    if (iter > max.iter) {
        result$convergence <- FALSE
    }
    class(result) <- "mle.wrappedcauchy"
    return(result)
} 

#############################################################
#                                                           #
#   print.mle.wrappednormal function                    #
#   Author: Claudio Agostinelli                         #
#   E-mail: claudio@unive.it                            #
#   Date: November, 19, 2003                                #
#   Version: 0.1-1                                        #
#                                                           #
#   Copyright (C) 2003 Claudio Agostinelli              #
#                                                           #
#############################################################

print.mle.wrappedcauchy <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    cat("mu: ")
    cat(format(x$mu, digits=digits), "\n")
    cat("\n")
    cat("rho: ")    
    cat(format(x$rho, digits=digits), "\n")
    if (!x$convergence) cat("\n The convergence is not achieved after the prescribed number of iterations \n")

    invisible(x)
}
