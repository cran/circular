intersect.modal.region <- function(x, ...) UseMethod("intersect.modal.region")

intersect.modal.region.default <- function(x, ...) .NotYetImplemented()

#############################################################
#
#	intersect.modal.region.circular
#       GNU General Public Licence 2.0 
#	Author: Claudio Agostinelli
#	E-mail: claudio.agostinelli@unitn.it
#	Date: November, 09, 2015
#	Version: 0.1
#
#	Copyright (C) 2015 Claudio Agostinelli
#
#############################################################

intersect.modal.region.circular <- function(x, breaks, z=NULL, q=0.95, bw, adjust = 1, type = c("K", "L"), kernel = c("vonmises", "wrappednormal"), na.rm = FALSE, step=0.01, eps.lower=10^(-4), eps.upper=10^(-4), ...) {
  breaks <- conversion.circular(breaks, units="radians", zero=0, rotation="counter", modulo="asis")    
  class(breaks) <- class(breaks)[class(breaks)!="circular"]
  attr(breaks, "circularp") <- NULL
  nrb <- nrow(breaks)
  mr <- modal.region.circular(x=x, z=z, q=q, bw=bw, adjust=adjust, type=type, kernel=kernel, na.rm=na.rm, step=step, eps.lower=eps.lower, eps.upper=eps.upper, ...)
  zeros <- mr$zeros
  nrz <- nrow(zeros)
  newbreaks <- areas <- list()
  tot <- 0
  for (i in 1:nrb) {
    newbreaks[[i]] <- matrix(0, nrow=0, ncol=2)  
    for (j in 1:nrz) {
      temp <- circular(IntersectIntervalRad(x=zeros[j,], y=breaks[i,]))
      newbreaks[[i]] <- rbind(newbreaks[[i]], temp)
    }
    if (nrow(newbreaks[[i]]) > 0) {
      areas[[i]] <- areas.region.circular(x=x, breaks=newbreaks[[i]], z=z, bw=bw, adjust=adjust, type=type, kernel=kernel, na.rm=na.rm, step=step, ...)
      tot <- tot + areas[[i]]$tot
    } else
      areas[[i]] <- NA
  }
  result <- list(tot=tot, areas=areas, breaks=newbreaks)
  return(result)
}

IntersectIntervalRad <- function(x, y) {
  x <- x%%(2*pi)
  y <- y%%(2*pi)
  if (x[1] <= x[2]) {
    if (y[1] <= y[2]) {
      if (x[2] < y[1] | y[2] < x[1])  
        res <- matrix(0, nrow=0, ncol=2)
      else
        res <- c(max(x[1], y[1]), min(x[2],y[2]))
    } else {
      res <- matrix(0, nrow=0, ncol=2)
      if (x[1] <= y[2]) {
        res <- rbind(res, c(x[1], min(y[2],x[2])))
      }
      if (x[2] >= y[1]) {
        res <- rbind(res, c(max(x[1],y[1]), x[2]))
      }
    }
  } else {
    if (y[1] <= y[2]) {
      res <- matrix(0, nrow=0, ncol=2)
      if (y[1] <= x[2]) {
        res <- rbind(res, c(y[1], min(x[2],y[2])))
      }
      if (y[2] >= x[1]) {
        res <- rbind(res, c(max(y[1],x[1]), y[2]))
      }
    } else {
      res <- rbind(c(0, min(x[2],y[2])),c(max(x[1],y[1]),2*pi))
      if (y[2] >= x[1])
        res <- rbind(res, c(x[1],y[2]))
      if (y[1] <= x[2])
        res <- rbind(res, c(y[1],x[2]))
    }
  }
  return(res)
}

if (FALSE) {
  library(gtools)  
  z <- pi/c(8,6,4,3)
  z <- permutations(n=4, r=4, v=z)
  for (i in 1:nrow(z)) {
    print(z[i,])
    print(IntersectIntervalRad(x=c(z[i,1],z[i,2]), y=c(z[i,3],z[i,4])))
  }
}

if (FALSE) {
  x <- rvonmises(100, circular(pi), 10)
  res <- intersect.modal.region(x, breaks=circular(matrix(c(pi,pi+pi/12), ncol=2)), bw=50)
  res$tot
  
  res <- intersect.modal.region(x, breaks=circular(matrix(c(pi,pi+pi/12, pi-pi/12, pi), ncol=2, byrow=TRUE)), bw=50)
  res$tot


  x <- rvonmises(100, circular(0), 10)
  res <- intersect.modal.region(x, breaks=circular(matrix(c(pi,pi+pi/12), ncol=2)), bw=50)
  res$tot
  
  res <- intersect.modal.region(x, breaks=circular(matrix(c(pi/12, 2*pi-pi/12), ncol=2, byrow=TRUE)), bw=50)
  res$tot

  res <- intersect.modal.region(x, breaks=circular(matrix(c(2*pi-pi/12, pi/12), ncol=2, byrow=TRUE)), bw=50)
  res$tot

  res <- intersect.modal.region(x, breaks=circular(matrix(c(2*pi-pi/12,2*pi), ncol=2, byrow=TRUE)), bw=50)
  res$tot

  res <- intersect.modal.region(x, breaks=circular(matrix(c(0,pi/12), ncol=2, byrow=TRUE)), bw=50)
  res$tot
}
