
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
#   Date: December, 5, 2005                                  #
#   Copyright (C) 2003 Claudio Agostinelli, David Andel     #
#   Copyright (C) 2005 Claudio Agostinelli     #
#                                                           #
#   Version 0.4                                           #
#############################################################

summary.circular <- function(object, ...) {
  if (is.matrix(object)) {
    return(summary.matrix(object, ...))
  }
  if (is.data.frame(object)) {
    return(summary.data.frame(object, ...))
  } else {
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
