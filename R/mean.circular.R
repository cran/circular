
###############################################################
#                                                             #
#       Original Splus: Ulric Lund                            #
#       E-mail: ulund@calpoly.edu                             #
#                                                             #
###############################################################

#############################################################
#                                                           #
#   mean.circular function                                  #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: December, 6, 2005                                   #
#   Version: 0.3-1                                            #
#                                                           #
#   Copyright (C) 2005 Claudio Agostinelli                  #
#                                                           #
#############################################################

mean.circular <- function(x, na.rm=FALSE, ...) {
   if (na.rm) 
       x <- x[!is.na(x)]
   if (length(x)==0) {
        warning("No observations (at least after removing missing values)")
        return(NULL)
   }
   x <- as.circular(x)
   xcircularp <- attr(x, "circularp")
   unitsp <- xcircularp$units

   if (any(is.na(x))) {
       circmean <- NA
   } else {
       if (unitsp=="degrees") 
           x <- x/180*pi

       sinr <- sum(sin(x))
       cosr <- sum(cos(x))

       if (sqrt((sinr^2 + cosr^2))/length(x) > .Machine$double.eps) {
           circmean <- atan2(sinr, cosr)
       } else {
           circmean <- NA
       }
   }
       
   if (unitsp=="degrees") 
       circmean <- circmean/pi*180

   attr(circmean, "circularp") <- xcircularp
   attr(circmean, "class") <- "circular"

   return(circmean)
}
