
#############################################################
#                                                           #
#       Original Splus: Ulric Lund                          #
#       E-mail: ulund@calpoly.edu                           #
#                                                           #
#############################################################

#############################################################
#                                                           #
#   rvonmises function                                      #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: August, 10, 2006                                  #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#   Version 0.2-4                                           #
#############################################################

rvonmises <- function(n, mu, kappa, control.circular=list()) {
   if (is.circular(mu)) {
      datacircularp <- circularp(mu)
   } else {
      datacircularp <- list(type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="counter")
   }
   dc <- control.circular
   if (is.null(dc$type))
      dc$type <- datacircularp$type
   if (is.null(dc$units))
      dc$units <- datacircularp$units
   if (is.null(dc$template))
      dc$template <- datacircularp$template
   if (is.null(dc$modulo))
      dc$modulo <- datacircularp$modulo
   if (is.null(dc$zero))
      dc$zero <- datacircularp$zero
   if (is.null(dc$rotation))
      dc$rotation <- datacircularp$rotation
   
   mu <- conversion.circular(mu, units="radians", zero=0, rotation="counter")
   attr(mu, "class") <- attr(mu, "circularp") <-  NULL  
   vm <- RvonmisesRad(n, mu, kappa)
   vm <- conversion.circular(circular(vm), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
   return(vm)
}

#RvonmisesRad <- function(n, mu, kappa) {
#   vm <- 1:n
#   a <- 1 + (1 + 4 * (kappa^2))^0.5
#   b <- (a - (2 * a)^0.5)/(2 * kappa)
#   r <- (1 + b^2)/(2 * b)
#  obs <- 1
#   while (obs <= n) {
#      U1 <- runif(1, 0, 1)
#      z <- cos(pi * U1)
#      f <- (1 + r * z)/(r + z)
#      c <- kappa * (r - f)
#      U2 <- runif(1, 0, 1)
#      if (c * (2 - c) - U2 > 0) {
#         U3 <- runif(1, 0, 1)
#         vm[obs] <- sign(U3 - 0.5) * acos(f) + mu
#         vm[obs] <- vm[obs] %% (2 * pi)
#         obs <- obs + 1
#      } else {
#         if (log(c/U2) + 1 - c >= 0) {
#            U3 <- runif(1, 0, 1)
#            vm[obs] <- sign(U3 - 0.5) * acos(f) + mu
#           vm[obs] <- vm[obs] %% (2 * pi)
#           obs <- obs + 1
#         }
#      }
#   }
#   return(vm)
#}

RvonmisesRad <- function(n, mu, kappa) {
   x <- vector(len = n)
   vm <- .C("rvm",
       as.double(x),
       as.integer(n),
       as.double(mu),
       as.double(kappa),
       PACKAGE="circular")[[1]] %% (2 * pi)
   return(vm)
}

#############################################################
#                                                           #
#   dvonmises function                                      #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: May, 10, 2006                                     #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#   Version 0.2                                             #
#############################################################

dvonmises <- function (x, mu, kappa) {
    x <- conversion.circular(x, units="radians", zero=0, rotation="counter")
    mu <- conversion.circular(mu, units="radians", zero=0, rotation="counter")

    attr(x, "class") <- attr(x, "circularp") <-  NULL
    attr(mu, "class") <- attr(mu, "circularp") <-  NULL    
  
    DvonmisesRad(x, mu, kappa)
}

DvonmisesRad <- function(x, mu, kappa) {
    return(1/(2 * pi * besselI(x = kappa, nu = 0, expon.scaled = TRUE)) * (exp(cos(x - mu) -1))^kappa)
}

#############################################################
#                                                           #
#   pvonmises function                                      #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: May, 10, 2006                                     #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#   Version 0.2                                             #
#############################################################

pvonmises <- function(q, mu, kappa, tol = 1e-020) {
     q <- conversion.circular(q, units="radians", zero=0, rotation="counter")
     mu <- conversion.circular(mu, units="radians", zero=0, rotation="counter")
     attr(q, "class") <- attr(q, "circularp") <-  NULL    
     attr(mu, "class") <- attr(mu, "circularp") <-  NULL

     PvonmisesRad(q, mu, kappa, tol)
}

PvonmisesRad <- function(q, mu, kappa, tol) {    
   q <- q %% (2 * pi)
   n <- length(q)
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
#   Date: May, 10, 2006                                     #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#   Version 0.2                                             #
#############################################################

dmixedvonmises <- function(x, mu1, mu2, kappa1, kappa2, p) {
    x <- conversion.circular(x, units="radians", zero=0, rotation="counter")
    mu1 <- conversion.circular(mu1, units="radians", zero=0, rotation="counter")
    mu2 <- conversion.circular(mu2, units="radians", zero=0, rotation="counter")
    
    attr(x, "class") <- attr(x, "circularp") <-  NULL
    attr(mu1, "class") <- attr(mu1, "circularp") <-  NULL
    attr(mu2, "class") <- attr(mu2, "circularp") <-  NULL
    
    DmixedvonmisesRad(x, mu1, mu2, kappa1, kappa2, p)
}

DmixedvonmisesRad <- function(x, mu1, mu2, kappa1, kappa2, p) {
   return(p/(2 * pi * besselI(x=kappa1, nu=0, expon.scaled = TRUE)) * (exp(cos(x - mu1) - 1))^kappa1 + (1 - p)/(2 * pi * besselI(x=kappa2, nu=0, expon.scaled = TRUE)) * (exp(cos(x - mu2) - 1))^kappa2)
}

#############################################################
#                                                           #
#   rmixedvonmises function                                 #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: August, 10, 2006                                  #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#   Version 0.2-4                                           #
#############################################################

rmixedvonmises <- function(n, mu1, mu2, kappa1, kappa2, p, control.circular=list()) {
   if (is.circular(mu1)) {
      datacircularp <- circularp(mu1)
   } else if  (is.circular(mu2)) {
      datacircularp <- circularp(mu2)
   } else {
      datacircularp <- list(type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="counter")
   }
   dc <- control.circular
   if (is.null(dc$type))
      dc$type <- datacircularp$type
   if (is.null(dc$units))
      dc$units <- datacircularp$units
   if (is.null(dc$template))
      dc$template <- datacircularp$template
   if (is.null(dc$modulo))
      dc$modulo <- datacircularp$modulo
   if (is.null(dc$zero))
      dc$zero <- datacircularp$zero
   if (is.null(dc$rotation))
      dc$rotation <- datacircularp$rotation
   
   mu1 <- conversion.circular(mu1, units="radians", zero=0, rotation="counter")
   mu2 <- conversion.circular(mu2, units="radians", zero=0, rotation="counter")
    
   attr(mu1, "class") <- attr(mu1, "circularp") <-  NULL
   attr(mu2, "class") <- attr(mu2, "circularp") <-  NULL

   vm <- RmixedvonmisesRad(n, mu1, mu2, kappa1, kappa2, p)
   vm <- conversion.circular(circular(vm), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)    
   return(vm)
}

RmixedvonmisesRad <- function(n, mu1, mu2, kappa1, kappa2, p) {
    result <- rep(NA, n)
    test <- runif(n)
    n1 <- sum(test < p)
    n2 <- n - n1
    res1 <- RvonmisesRad(n1, mu1, kappa1)
    res2 <- RvonmisesRad(n2, mu2, kappa2)
    result[test < p] <- res1
    result[test >= p] <- res2
    return(result)
}
