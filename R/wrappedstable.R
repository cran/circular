#############################################################
#                                                           #
#   rwrappedstable function                                 #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: May, 29, 2006                                     #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#   Version 0.2-2                                           #
#############################################################

rwrappedstable <- function(n,  scale=1, index, skewness, control.circular=list()) {
    dc <- control.list
    result <- rstable(n=n, scale=scale, index=index, skewness=skewness) %% (2 * pi)
    result <- conversion.circular(circular(result), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
    return(result)
}

