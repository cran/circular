
###############################################################
#                                                             #
#       Original Splus: Ulric Lund                            #
#       E-mail: ulund@calpoly.edu                             #
#                                                             #
###############################################################

#############################################################
#                                                           #
#   change.point function                                   #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: July, 24, 2003                                    #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1                                             #
#############################################################

change.point <- function(x) {
        x <- as.circular(x)
        xcircularp <- circularp(x)
        units <- xcircularp$units
        x <- conversion.circular(x, units="radians")
    phi <- function(x) {
        arg <- A1inv(x)
        if(besselI(x=arg, nu=0, expon.scaled = FALSE) != Inf)
            result <- x * A1inv(x) - log(besselI(x=arg, nu=0, expon.scaled = FALSE))
        else result <- x * A1inv(x) - (arg + log(1/sqrt(2 * pi * arg) * (1 + 1/(8 * arg) + 9/(128 * arg^2) + 225/(1024 * arg^3))))
        result
    }
    n <- length(x)
    rho <- rho.circular(x)
    R1 <- c(1:n)
    R2 <- c(1:n)
    V <- c(1:n)
    for(k in 1:(n - 1)) {
        R1[k] <- rho.circular(x[1:k]) * k
        R2[k] <- rho.circular(x[(k + 1):n]) * (n - k)
        if(k >= 2 & k <= (n - 2)) {
            V[k] <- k/n * phi(R1[k]/k) + (n - k)/n * phi(R2[k]/(n - k))
        }
    }
    R1[n] <- rho * n
    R2[n] <- 0
    R.diff <- R1 + R2 - rho * n
    rmax <- max(R.diff)
    rave <- mean(R.diff)
    k.r <- (1:n)[R.diff == max(R.diff)]
    V <- V[2:(n - 2)]
    if(n > 3) {
        tmax <- max(V)
        tave <- mean(V)
        k.t <- (1:(n - 3))[V == max(V)] + 1
    }
    else stop("Sample size must be at least 4")
    return(list(n=n, rho=rho, rmax=rmax, k.r=k.r, rave=rave, tmax=tmax, k.t=k.t, tave=tave))
}
