###############################################################
#                                                             #
#       Original Splus: Ulric Lund                            #
#       E-mail: ulund@calpoly.edu                             #
#                                                             #
###############################################################


#############################################################
#                                                           #
#   range.circular function                                 #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: December, 23, 2003                                #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1-1                                           #
#############################################################

range.circular <- function(x, test = FALSE, na.rm=FALSE, finite=FALSE, ...) {
       x <- as.circular(x)
       units <- circularp(x)$units
       x <- conversion.circular(x, units="radians")
       if (finite) 
           x <- x[is.finite(x)]
        else if (na.rm) 
           x <- x[!is.na(x)]
        if (length(x)) 
            c(min(x), max(x))
        else c(NA, NA)
    x <- sort(x %% (2*pi))
    n <- length(x)
    spacings <- c(diff(x), x[1] - x[n] + 2*pi)
    range <- 2*pi - max(spacings)
        if (units=="degrees") 
            rangenew <- range/pi*180
        else
            rangenew <- range
        result <- rangenew
    if(test == TRUE) {
        stop <- floor(1/(1 - range/(2*pi)))
        index <- c(1:stop)
        sequence <- ((-1)^(index - 1)) * exp(log(gamma(n + 1)) - log(gamma(index + 1)) - log(gamma(n - index + 1))) * (1 - index * (1 - range/(2 * pi)))^(n - 1)
        p.value <- sum(sequence)
        result <- list(range=rangenew, p.value=p.value)
    }
    result
}
