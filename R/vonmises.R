
###############################################################
#                                                             #
#       Original Splus: Ulric Lund                            #
#       E-mail: ulund@calpoly.edu                             #
#                                                             #
###############################################################

#############################################################
#                                                           #
#   rvonmises function                                      #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: July, 21, 2003                                    #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1                                             #
#############################################################

rvonmises <- function(n, mu, kappa, units=c("radians", "degrees"), ...) {
    units <- match.arg(units)
    if (units=="degrees") {
        mu <- mu/180*pi
    }

    vm <- 1:n
    a <- 1 + (1 + 4 * (kappa^2))^0.5
    b <- (a - (2 * a)^0.5)/(2 * kappa)
    r <- (1 + b^2)/(2 * b)
    obs <- 1
    while (obs <= n) {
       U1 <- runif(1, 0, 1)
       z <- cos(pi * U1)
       f <- (1 + r * z)/(r + z)
       c <- kappa * (r - f)
       U2 <- runif(1, 0, 1)
       if (c * (2 - c) - U2 > 0) {
           U3 <- runif(1, 0, 1)
           vm[obs] <- sign(U3 - 0.5) * acos(f) + mu
           vm[obs] <- vm[obs] %% (2 * pi)
           obs <- obs + 1
       } else {
           if (log(c/U2) + 1 - c >= 0) {
           U3 <- runif(1, 0, 1)
           vm[obs] <- sign(U3 - 0.5) * acos(f) + mu
           vm[obs] <- vm[obs] %% (2 * pi)
           obs <- obs + 1
           }
       }
    }
    if (units=="degrees") vm <- vm/pi*180
    vm <- circular(vm, units=units, ...)
    return(vm)
}

#############################################################
#                                                           #
#   dvonmises function                                      #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: July, 21, 2003                                    #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1                                             #
#############################################################

dvonmises <- function (x, mu, kappa) {
    x <- as.circular(x)
    xcircularp <- circularp(x)
    units <- xcircularp$units
    x <- conversion.circular(x, units="radians")
    n <- length(x)
    attr(x, "class") <- attr(x, "circularp") <-  NULL
    
    if (units=="degrees") {
        mu <- mu/180*pi
    }  
  
    return(1/(2 * pi * besselI(x = kappa, nu = 0, expon.scaled = TRUE)) * (exp(cos(x - mu) -1))^kappa)
}

#############################################################
#                                                           #
#   pvonmises function                                      #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: July, 24, 2003                                    #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1-1                                           #
#############################################################

pvonmises <- function(q, mu, kappa, tol = 1e-020) {

     q <- as.circular(q)
     qcircularp <- circularp(q)
     units <- qcircularp$units
     q <- conversion.circular(q, units="radians")
     n <- length(q)
     attr(q, "class") <- attr(q, "circularp") <-  NULL
    
     if (units=="degrees") {
         mu <- mu/180*pi
     }  
     attr(mu, "class") <- attr(mu, "circularp") <-  NULL
 
     q <- q %% (2 * pi)
     mu <- mu %% (2 * pi)
     pvm.mu0 <- function(q, kappa, tol) {
        flag <- TRUE
        p <- 1
        sum <- 0
        while (flag) {
           term <- (besselI(x=kappa, nu=p, expon.scaled = FALSE) * sin(p * q))/p
           sum <- sum + term
           p <- p + 1
           if (abs(term) < tol)
           flag <- FALSE
    }
    return(q/(2 * pi) + sum/(pi * besselI(x=kappa, nu=0, expon.scaled = FALSE)))
     }
     result <- rep(NA, n)
     if (mu == 0) {
         for (i in 1:n) {
          result[i] <- pvm.mu0(q[i], kappa, tol)
         }
     } else {
         for (i in 1:n) {
           
         if (q[i] <= mu) {
             upper <- (q[i] - mu) %% (2 * pi)
             if (upper == 0)
             upper <- 2 * pi
             lower <- ( - mu) %% (2 * pi)
             result[i] <- pvm.mu0(upper, kappa, tol) - pvm.mu0(lower, kappa, tol)
             } else {
             upper <- q[i] - mu
             lower <- mu %% (2 * pi)
             result[i] <- pvm.mu0(upper, kappa, tol) + pvm.mu0(lower, kappa, tol)
            }
         }     
     }
    return(result)
}

#############################################################
#                                                           #
#   dmixedvonmises function                                 #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: July, 21, 2003                                    #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1                                             #
#############################################################

dmixedvonmises <- function(x, mu1, mu2, kappa1, kappa2, p) {
    x <- as.circular(x)
    xcircularp <- circularp(x)
    units <- xcircularp$units
    x <- conversion.circular(x, units="radians")
    n <- length(x)
    attr(x, "class") <- attr(x, "circularp") <-  NULL
    
    if (units=="degrees") {
        mu1 <- mu1/180*pi
        mu2 <- mu2/180*pi
    }  
  
    return(p/(2 * pi * besselI(x=kappa1, nu=0, expon.scaled = TRUE)) * (exp(cos(x - mu1) - 1))^kappa1 + (1 - p)/(2 * pi * besselI(x=kappa2, nu=0, expon.scaled = TRUE)) * (exp(cos(x - mu2) - 1))^kappa2)
}

#############################################################
#                                                           #
#   rmixedvonmises function                                 #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: July, 21, 2003                                    #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1                                             #
#############################################################

rmixedvonmises <- function(n, mu1, mu2, kappa1, kappa2, p, units=c("radians", "degrees"), ...) {
    units <- match.arg(units)
    result <- rep(NA, n)
    test <- runif(n)
    n1 <- sum(test < p)
    n2 <- n - n1
    res1 <- rvonmises(n1, mu1, kappa1, units=units, ...)
    res2 <- rvonmises(n2, mu2, kappa2, units=units, ...)
    result[test < p] <- res1
    result[test >= p] <- res2
    return(result)
}
