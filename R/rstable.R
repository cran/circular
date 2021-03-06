###############################################################
#       rstable function                                      #
#       Date: January, 22, 2002                               #
#       Version: 0.1                                          #
#                                                             #
###############################################################
#                                                             #
#   This  R code is based on C functions gsl_ran_levy and     #
#     gsl_ran_levy_skew from GNU Scientifi Library            #
#      copyrighted under GNU general license by               #
#     James Theiler, Brian Gough and Keith Briggs.            #
#                                                             #
###############################################################     
#   Here the original comments in the code:
#
#   The stable Levy probability distributions have the form
#
#   p(x) dx = (1/(2 pi)) \int dt exp(- it x - |c t|^alpha)
#
#   with 0 < alpha <= 2. 
#
#   For alpha = 1, we get the Cauchy distribution
#   For alpha = 2, we get the Gaussian distribution with sigma = sqrt(2) c.
#
#   Fromn Chapter 5 of Bratley, Fox and Schrage "A Guide to
#   Simulation". The original reference given there is,
#
#   J.M. Chambers, C.L. Mallows and B. W. Stuck. "A method for
#   simulating stable random variates". Journal of the American
#   Statistical Association, JASA 71 340-344 (1976).
#
#   The following routine for the skew-symmetric case was provided by
#   Keith Briggs.
#
#   The stable Levy probability distributions have the form
#
#   2*pi* p(x) dx
#
#     = \int dt exp(mu*i*t-|sigma*t|^alpha*(1-i*beta*sign(t)*tan(pi*alpha/2))) for alpha!=1
#     = \int dt exp(mu*i*t-|sigma*t|^alpha*(1+i*beta*sign(t)*2/pi*log(|t|)))   for alpha==1
#
#   with 0<alpha<=2, -1<=beta<=1, sigma>0.
#
#   For beta=0, sigma=c, mu=0, we get gsl_ran_levy above.
#
#   For alpha = 1, beta=0, we get the Lorentz distribution
#   For alpha = 2, beta=0, we get the Gaussian distribution
#
#   See A. Weron and R. Weron: Computer simulation of Levy alpha-stable 
#   variables and processes, preprint Technical University of Wroclaw.
#   http://www.im.pwr.wroc.pl/~hugo/Publications.html
#
###############################################################################

rstable <- function(n, scale = 1, index = stop("no index arg"), skewness = 0) {
##    alpha <- index
##    beta <- skewness
    if (index > 2 | index <= 0)
       stop("rstable is not define for index outside the interval 0 < index <= 2\n")
    if (skewness > 1 | skewness < -1)
       {stop("rstable is not define for skewness outside the interval -1 <= skewness <= 1\n")}
    if (skewness==0) {
## cauchy case
       if (index == 1) { 
          return(scale*rcauchy(n, location = 0, scale = 1))
       }
## gaussian case
       if (index == 2) { 
          return(rnorm(n, mean = 0, sd = sqrt(2)*scale)) 
       }
## general case
       rngstab <- vector(length=0)
       for (i in 1:n) {
          u <- 0 
          while (u == 0 | u == 1) {
              u <-  pi * (runif(1, min=0, max=1) - 0.5)
          }
          v <- 0
          while (v == 0) {
              v <- rexp(1,rate=1)   
          }
          t <-  sin (index * u) / (cos (u)^(1 / index))
          s <- (cos ((1 - index) * u) / v)^((1 - index) / index)
          rngstab <- c(rngstab, t*s)
       }
       return (scale * rngstab);
    } else {
       rngstab <- vector(length=0)
       for (i in 1:n) {
          u <- 0 
          while (u == 0 | u == 1) {
              u <-  pi * (runif(1, min=0, max=1) - 0.5)
          }
          v <- 0
          while (v == 0) {
              v <- rexp(1,rate=1)   
          }
          if (index == 1) {
             X <-  (((pi/2) + skewness * u) * tan (u) - skewness * log ((pi/2) * v * cos (u) / ((pi/2) + skewness * u))) / (pi/2)
             rngstab <- c(rngstab, (scale * (X + skewness * log (scale) / (pi/2))))
           } else {
              t <- skewness * tan ((pi/2) * index)
              B <- atan(t) / index
              S <-  (1 + t * t)^(1/(2 * index))
              X <-  S * sin (index * (u + B)) / (cos (u)^(1 / index)) * (cos (u - index * (u + B)) / v)^((1 - index) / index)
              rngstab <- c(rngstab, (scale * X))
           }
        }
        return(rngstab)
    }
}

