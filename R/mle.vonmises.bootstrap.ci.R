
###############################################################
#                                                             #
#       Original Splus: Ulric Lund                            #
#       E-mail: ulund@calpoly.edu                             #
#                                                             #
###############################################################

#############################################################
#                                                           #
#   mle.vonmises.bootstrap.ci function                      #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: April, 11, 2005                                   #
#   Copyright (C) 2005 Claudio Agostinelli                  #
#                                                           #
#   Version 0.2                                             #
#############################################################

mle.vonmises.bootstrap.ci <- function(x, mu, bias = FALSE, alpha = 0.05, reps = 1000) {
  
  # Handling missing values
  x <- na.omit(x)
  if (length(x)==0) {
      warning("No observations (at least after removing missing values)")
      return(NULL)
  }
  
  if (require(boot)) {    
      x <- as.circular(x)
      xcircularp <- circularp(x)
      units <- xcircularp$units
      x <- conversion.circular(x, units="radians")

      if (missing(mu)) {
          sinr <- sum(sin(x))
          cosr <- sum(cos(x))
          mu <- atan(sinr, cosr)
      } else {
          attr(mu, "circularp") <- xcircularp
          attr(mu, "class") <- "circular"
          mu <- conversion.circular(mu, units="radians")
      }
    
      mle.vonmises.mu <- function(x, i) {
          sinr <- sum(sin(x[i]))
          cosr <- sum(cos(x[i]))
          mu <- atan(sinr, cosr)
          return(mu)
      }

      mle.vonmises.kappa <- function(x, i, mu, bias) {
          n <- length(x[i])
          V <- mean(cos(x[i] - mu))
          if (V > 0) {
              kappa <- A1inv(V)
          } else {
              kappa <- 0
          }
          if (bias == TRUE) {
              if (kappa < 2) {
                  kappa <- max(kappa - 2 * (n * kappa)^-1, 0)
              } else {
                  kappa <- ((n - 1)^3 * kappa)/(n^3 + n)
              }
          }
          return(kappa)
      }
      
      mean.bs <- boot(data = x, statistic = mle.vonmises.mu, R = reps, stype="i")

      mean.reps <- mean.bs$t
      mean.reps <- sort(mean.reps %% (2 * pi))
      B <- reps
      spacings <- c(diff(mean.reps), mean.reps[1] - mean.reps[B] + 2 * pi)
      max.spacing <- (1:B)[spacings == max(spacings)]
      off.set <- 2 * pi - mean.reps[max.spacing + 1]
      if (max.spacing != B)
      mean.reps2 <- mean.reps + off.set
      else mean.reps2 <- mean.reps
      mean.reps2 <- sort(mean.reps2 %% (2 * pi))
      mean.ci <- quantile(mean.reps2, c(alpha/2, 1 - alpha/2))
      if (max.spacing != B)
      mean.ci <- mean.ci - off.set
    
      
      kappa.bs <- boot(data = x, statistic = mle.vonmises.kappa, R = reps, stype="i", mu=mu, bias = bias)
      kappa.reps <- kappa.bs$t

      kappa.ci <- quantile(kappa.reps, c(alpha/2, 1 - alpha/2))

      if (units=="degrees") {
          mean.reps <- mean.reps/pi*180
          mean.ci <- mean.ci/pi*180
      }
      attr(mean.reps, "circularp") <- xcircularp
      attr(mean.reps, "class") <- "circular"
      attr(mean.ci, "circularp") <- xcircularp
      attr(mean.ci, "class") <- "circular"

      result <- list()
      result$call <- match.call()
      result$mu.ci <- mean.ci
      result$mu <- c(mean.reps)
      result$kappa.ci <- kappa.ci
      result$kappa <- c(kappa.reps)
      result$alpha <- alpha
      class(result) <- "mle.vonmises.bootstrap.ci"
      return(result)

   } else {
       stop("To use this function you have to install the package 'boot' \n")
   }

}


#############################################################
#                                                           #
#   print.mle.vonmises.bootstrap.ci function                #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: September, 17, 2003                               #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1                                             #
#############################################################

print.mle.vonmises.bootstrap.ci <- function(x, ...) {
    cat("Bootstrap Confidence Intervals for Mean Direction and Concentration", "\n")
    cat("Confidence Level:  ", round(100 * (1 - x$alpha),2), "%", "\n")
    cat("Mean Direction:           ", "Low =", round(x$mu.ci[1], 2), "  High =", round(x$mu.ci[2], 2), "\n")
    cat("Concentration Parameter:  ", "Low =", round(x$kappa.ci[1], 2), "  High =", round(x$kappa.ci[2], 2), "\n")

}
        
