#############################################################
#                                                           #
#   rpnorm function                                         #
#   Author: Claudio Agostinelli                             #
#   Email: claudio.agostinelli@unitn.it                     #
#   Date: December, 08, 2016                                #
#   Copyright (C) 2016 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1                                             #
#############################################################

rpnorm <- function(n, mu, sigma, control.circular=list()) {
  if (missing(mu) || length(mu)!=2)
    stop("the mean direction parameter 'mu' is mandatory and it must have length 1")
  if (missing(sigma) || dim(sigma)!=c(2,2))
    stop("the variance matrix parameter 'Sigma' is mandatory and it must be a matrix of dimension 2 by 2")   
  datacircularp <- list(type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="counter")
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
  mu <- as.vector(mu)
  pn <- RpnormRad(n, mu, sigma)
  pn <- conversion.circular(circular(pn), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
  return(pn)
}

RpnormRad <- function(n, mu, sigma) {
  x <- rmvnorm(n=n, mean=mu, sigma=sigma)
  x <- apply(x, 1, function(z) z/sqrt(z[1]^2+z[2]^2))
  theta <- apply(x, 2, function(z) atan2(z[2],z[1]))
  theta <- theta%%(2*pi)
  return(theta)
}

#############################################################
#                                                           #
#   dpnorm function                                         #
#   Author: Claudio Agostinelli                             #
#   Email: claudio.agostinelli@unitn.it                     #
#   Date: December, 08, 2016                                #
#   Copyright (C) 2016 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1                                             #
#############################################################

dpnorm <- function (x, mu, sigma, log=FALSE) {
  if (missing(mu) || length(mu)!=2)
    stop("the mean direction parameter 'mu' is mandatory and it must have length 2")
  if (missing(sigma) || dim(sigma)!=c(2,2))
    stop("the variance matrix parameter 'sigma' is mandatory and it must be a matrix of dimension 2 by 2")   
  if (!is.logical(log))
    stop("'log' must be logical")
  x <- conversion.circular(x, units="radians", zero=0, rotation="counter")
  mu <- as.vector(mu)
  attr(x, "class") <- attr(x, "circularp") <-  NULL
  DpnormRad(x, mu, sigma, log)
}

DpnormRad <- function(x, mu, sigma, log=FALSE) {
  ### Wang and Gelfand (2013) Statistical Methodology  
  sigma1 <- sqrt(sigma[1,1])
  sigma2 <- sqrt(sigma[2,2])
  rho <- sigma[1,2]/(sigma1*sigma2)
  a <- 1/(sigma1*sigma2*sqrt(1-rho^2))
  C <- a^2*(sigma2^2*cos(x)^2-rho*sigma1*sigma2*sin(2*x)+sigma1^2*sin(x)^2)
  D <- a^2*(mu[1]*sigma2*(sigma2*cos(x)-rho*sigma1*sin(x))+mu[2]*sigma1*(sigma1*sin(x)-rho*sigma2*cos(x)))/sqrt(C)
  den <- dmvnorm(mu, c(0,0), sigma)/C+a*D*pnorm(D)*dnorm(a*(mu[1]*sin(x)-mu[2]*cos(x))/sqrt(C))/C
  if (log) {
    den <- log(den)  
  }
  return(den)
}
