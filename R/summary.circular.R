
###############################################################
#                                                             #
#       Original Splus: Ulric Lund                            #
#       E-mail: ulund@calpoly.edu                             #
#                                                             #
###############################################################

#############################################################
#                                                           #
#   summary.circular                                        #
#   Authors: Claudio Agostinelli, David Andel               #
#   Email: claudio@unive.it, andel@ifi.unizh.ch             #
#   Date: August, 01, 2003                                  #
#   Copyright (C) 2003 Claudio Agostinelli, David Andel     #
#                                                           #
#   Version 0.3-1                                           #
#############################################################

summary.circular <- function(object, ...) {
  if (is.matrix(object)) {
    return(summary.matrix(object, ...))
  }
  else {
    nas <- is.na(object)
    object <- object[!nas]
    n <- length(object)
    result <- c(n, mean.circular(object), rho.circular(object))
    names(result) <- c("n", "Mean", "Rho")
    if(any(nas))
      c(result, "NA's" = sum(nas))
    else result
  }
}
