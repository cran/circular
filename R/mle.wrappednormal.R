#############################################################
#                                                           #
#   mle.wrappednormal function                              #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: July, 31, 2003                                    #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1-2                                           #
#############################################################

mle.wrappednormal <- function(x, mu, rho, sd, K, tol=1e-5, min.sd=1e-3, min.k=10, max.iter=100, verbose=FALSE) {

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
    }
    est.rho <- FALSE
    if (missing(sd)) {
        if (missing(rho)) {
            sd <- sqrt(-2*log(sqrt(sinr^2 + cosr^2)/n))
            if (is.na(sd) || sd < min.sd) sd <- min.sd
            est.rho <- TRUE
        } else {
            sd <- sqrt(-2*log(rho))
        }
    }
     
    xdiff <- 1+tol
    iter <- 0
    if (missing(K)) {
        range <- max(mu, x) - min(mu, x)
        K <- (range+6*sd)%/%(2*pi)+1
        K <- max(min.k, K)
    }
    
    while (xdiff > tol & iter <= max.iter) {
           iter <- iter + 1
           mu.old <- mu
           sd.old <- sd

           z <- .Fortran("mlewrpno",
                    as.double(x),
                    as.double(mu),
                    as.double(sd),
                    as.integer(n),
                    as.integer(K),
                    as.integer(est.mu),
                    as.integer(est.rho),
                    w=double(n),
                    wk=double(n),
                    wm=double(n),
                    PACKAGE="circular"
           )
           w <- z$w
           wk <- z$wk
           wm <- z$wm
           
           if (est.mu) {
               mu <- sum(x)/n
               if (any(wk!=0)) {
                   mu <- mu + 2*pi*mean(wk[wk!=0]/w[wk!=0])
               }
           }
           if (est.rho) {
               if  (any(wm!=0)) {
                    sd <- sqrt(sum(wm[wm!=0]/w[wm!=0])/n)
               } else {
                    sd <- min.sd
               }
           }

           if (verbose) {
               cat("mu: ", mu, "\n")
               cat("rho: ", exp(-sd^2/2), "\n")              
               cat("sd: ", sd, "\n")
           }
           xdiff <- max(abs(mu - mu.old), abs(sd - sd.old))
    }

    rho <- exp(-sd^2/2)

    if (units=="degrees") {
        mu <- mu/pi*180
        sd <- sd/pi*180
    }
    
    attr(mu, "circularp") <- xcircularp
    attr(mu, "class") <- "circular"
    
    result <- list()
    
    result$call <- match.call()
    result$mu <- mu
    result$rho <- rho
    result$sd <- sd
    result$est.mu <- est.mu
    result$est.rho <- est.rho
    result$convergence <- TRUE
    if (iter > max.iter) {
        result$convergence <- FALSE
    }
    class(result) <- "mle.wrappednormal"
    return(result)
}

#############################################################
#                                                           #
#	print.mle.wrappednormal function                    #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: November, 19, 2003                                #
#	Version: 0.1-2                                      #
#                                                           #
#	Copyright (C) 2003 Claudio Agostinelli              #
#                                                           #
#############################################################

print.mle.wrappednormal <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    cat("mu: ")
    cat(format(x$mu, digits=digits), "\n")
    cat("\n")
    cat("rho: ")    
    cat(format(x$rho, digits=digits), "\n")
    cat("\n")
    cat("sd: ")       
    cat(format(x$sd, digits=digits), "\n")
    cat("\n")   
    if (!x$est.mu) cat("mu is known\n")
    if (!x$est.rho) {
        cat("rho and sd are known\n")
    }
    if (!x$convergence) cat("\nThe convergence is not achieved after the prescribed number of iterations \n")
    invisible(x)
}

