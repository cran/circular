#############################################################
#
#	totalvariation.circular
#       GNU General Public Licence 2.0 
#	Author: Claudio Agostinelli
#	E-mail: claudio@unive.it
#	Date: May, 4, 2011
#	Version: 0.4
#
#	Copyright (C) 2011 Claudio Agostinelli
#
#############################################################

totalvariation.circular <- function(x, y, z=NULL, q=0.95, bw, adjust = 1, type = c("K", "L"), kernel = c("vonmises", "wrappednormal"), na.rm = FALSE, step=0.001, eps.lower=10^(-4), eps.upper=10^(-4), ...) {
  if (is.null(z))
    z <- circular(seq(0,2*pi+step,step))
  if (!is.circular(x)) {
    x <- circular(x)
    cat("'x' is coerced to circular object assuming default values for the 'circular' function\n")
  }
  if (!is.circular(y)) {
    y <- circular(y)
    cat("'y' is coerced to circular object assuming default values for the 'circular' function\n")    
  }    
  modalx <- modal.region.circular(x=x, z=z, q=q, bw=bw, adjust=adjust, type=type, kernel=kernel, na.rm=na.rm, step=step, eps.lower=eps.lower, eps.upper=eps.upper, ...)
  modaly <- modal.region.circular(x=y, z=z, q=q, bw=bw, adjust=adjust, type=type, kernel=kernel, na.rm=na.rm, step=step, eps.lower=eps.lower, eps.upper=eps.upper, ...)
  zerosx <- modalx$zeros
  areasx <- modalx$areas$tot
  zerosy <- modaly$zeros
  nx <- nrow(zerosx)
  ny <- nrow(zerosy)
  areasy <- modaly$areas$tot
  denx <- modalx$density
  deny <- modaly$density  
  denx <- denx$y/areasx
  deny <- deny$y/areasy
  denx[z < zerosx[1,1]] <- 0
  if (nx > 1) {
    for (i in 1:(nx-1)) {
      denx[z > zerosx[i,2] & z < zerosx[i+1,1]] <- 0
    }
  }
  denx[z > zerosx[nx, 2]] <- 0
  deny[z < zerosy[1,1]] <- 0
  if (ny > 1) {
    for (i in 1:(ny-1)) {
      deny[z > zerosy[i,2] & z < zerosy[i+1,1]] <- 0
    }
  }
  deny[z > zerosy[ny, 2]] <- 0  
  den <- ifelse(denx > deny, denx-deny, 0) 
  denmax <- approxfun(x=z, y=den)
  tv <- 0
##  byhand <- 0    
  for (i in 1:nx) {
##    byhand <- byhand + step*sum(den[z >= zerosx[i,1] & z <= zerosx[i,2]])
    tv <- tv +  integrate(denmax, lower=zerosx[i,1], upper=zerosx[i,2])$value
  }
  result <- list()
  result$tv <- tv
##  result$byhand <- byhand
  result$q <- q
  result$x <- modalx
  result$y <- modaly
  class(result) <- 'totalvariation.circular'
  return(result)
}
