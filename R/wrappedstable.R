#############################################################
#                                                           #
#   rwrappedstable function                                 #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: September, 22, 2003                               #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#   Version 0.2                                             #
#############################################################

rwrappedstable <- function(n,  scale=1, index, skewness, units=c("radians", "degrees"), ...) {
    units <- match.arg(units)
#    if (units=="degrees") {
#        scale <- scale/180*pi
#    }
    result <- rstable(n=n, scale=scale, index=index, skewness=skewness) %% (2 * pi)
    if (units=="degrees") result <- result/pi*180
    result <- circular(result, units=units, ...)
    return(result)
}

