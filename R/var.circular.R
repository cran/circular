
var <- function(x, ...) UseMethod("var")

var.default <- function(x, y = NULL, na.rm = FALSE, use, ...) base::var(x=x, y=y, na.rm=na.rm, use=use)

#var.matrix <- function(x, ...) {
#    apply(x, 2, var, ...)
#}

var.data.frame <- function(x, ...) {
    sapply(x, var, ...)
}


###############################################################
#                                                             #
#       Original Splus: Ulric Lund                            #
#       E-mail: ulund@calpoly.edu                             #
#                                                             #
###############################################################

#############################################################
#                                                           #
#   var.circular function                                   #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: November, 19, 2003                                #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1                                             #
#############################################################

var.circular <- function (x, na.rm = FALSE, only.var=TRUE, ...)  {
    x <- as.circular(x)
    x <- conversion.circular(x, units="radians")

    if (is.matrix(x)) {
        apply(x, 2, var.circular, na.rm=na.rm, only.var=only.var)
    } else {

       if (na.rm) 
           x <- x[!is.na(x)]
       if (any(is.na(x))) return(NA)

       n <- length(x)
       c <- sum(cos(x))
       s <- sum(sin(x))
       r <- sqrt(c^2 + s^2)
       rbar <- r/n
       circvar <- 1 - rbar
       if (only.var) {
           return(circvar)
       } else {
           return(c(n,r,rbar,circvar))
       }
   }
}