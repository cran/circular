#############################################################
#                                                           #
#   axis.circular function                                  #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: October, 07, 2003                                 #
#   Version: 0.3-2                                          #
#                                                           #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################
 
axis.circular <- function(at, labels,  units = c("radians", "degrees"), template=c("none", "geographics"), modulo = c("asis", "2pi", "pi"), zero=0, rotation = c("counter", "clock"), tick=TRUE, lty, lwd, cex, col, font, tcl=0.025, tcl.text=0.125, digits=2) {

  units <- match.arg(units)
  template <- match.arg(template)
  modulo <- match.arg(modulo)
  rotation <- match.arg(rotation)

  if (missing(at)) {
      if (units=="radians") {
          at <- c(0, pi/2, pi, 3/2*pi)
      } else {
          at <- c(0, 90, 180, 270)
      }   
  }

  at <- as.circular(at, units=units, template=template, modulo=modulo, zero=zero, rotation=rotation)
  atcircularp <- circularp(at)
  zero <- atcircularp$zero
  rotation <- atcircularp$rotation
  atasis <- at
  attr(atasis, "circularp") <- attr(atasis, "class") <- NULL
  at <- conversion.circular(at, units="radians")
  attr(at, "circularp") <- attr(at, "class") <- NULL

  attext <- round(at/pi, digits=digits)
  
  if (missing(cex)) cex <- par("cex.axis")
  if (missing(col)) col <- par("col.axis")
  if (missing(font)) font <- par("font.axis")
  if (missing(lty)) lty <- par("lty")
  if (missing(lwd)) lwd <- par("lwd")

  if (missing(labels) & all(at==c(0, pi/2, pi, 3/2*pi))) { 
      if (template=="geographics") {
          labels <- c("N", "E", "S", "W")         
      } else {
          if (units=="radians") {
              labels <- c("0", expression(frac(pi,2)), expression(pi), expression(frac(3*pi,2)))
          } else {
              labels <- c("0", "90", "180", "270")      
          }
      }   
  } else {
      if (missing(labels)) {
          if (units!="radians") {
              labels <- as.character(round(atasis, digits=digits))
          }
      }
  }

if (!missing(labels) && length(at)!=length(labels)) stop("'at' and 'labels' must have the same length")
  
  if (rotation=="clock") at <- -at 
  at <- at + zero
        
  r <- 1+tcl*c(-1/2,1/2)
  r.l <- 1-tcl.text 
  z <- cos(at)
  y <- sin(at)
  
  for (i in 1:length(at)) {
       if (tick) {
           lines.default(z[i]*r, y[i]*r, col=col, lty=lty, lwd=lwd)
       }
       if (missing(labels) & units=="radians") {
              labeltext <- substitute(at*pi, list(at=attext[i]))
       } else {
              labeltext <- labels[i]
       }
       text.default(z[i]*r.l, y[i]*r.l, labeltext, cex=cex, col=col)    
  }
     
  text(0, 0, "+", cex=1)
}

#############################################################
#                                                           #
#   ticks.circular function                                 #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: September, 22, 2003                               #
#   Version: 0.2-1                                          #
#                                                           #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################
 
ticks.circular <- function(x, template=c("none", "geographics"), zero, rotation, tcl=0.025, col, ...) {

  template <- match.arg(template)
  if (missing(col)) col <- par("col")

  x <- as.circular(x, template=template)
  xcircularp <- attr(x, "circularp")
  xzero <- xcircularp$zero
  xrotation <- xcircularp$rotation
  attr(x, "circularp") <- attr(x, "class") <- NULL

  if (missing(zero)) {
      if (template=="geographics") {
          zero <- pi/2
      } else {
          zero <- xzero
      }
  }
  
  if (missing(rotation)) {
      if (template=="geographics") {
          rotation <- "clock"
      } else {
          rotation <- xrotation
      }
  }
  
  if (rotation=="clock") x <- -x 
  x <- x + zero
        
  r <- 1+tcl*c(-1/2,1/2)
  z <- cos(x)
  y <- sin(x)

  for (i in 1:length(x)) {
       lines.default(z[i]*r, y[i]*r, col=col, ...)
  }
}

#############################################################
#                                                           #
#   plot.circular function                                  #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: November, 18, 2003                                #
#   Version: 0.2-5                                          #
#                                                           #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################
 
plot.circular <- function(x, pch = 16, cex = 1, stack = FALSE, axes = TRUE, sep = 0.025, shrink = 1, bins, ticks = FALSE, tcl=0.025, tcl.text=0.125, col, tol = 0.04, uin, xlim=c(-1, 1), ylim=c(-1, 1), main=NULL, digits=2, ...) {

    if (is.matrix(x) | is.data.frame(x)) {
        nseries <- ncol(x)
    } else {
        nseries <- 1
    }
    xx <- as.data.frame(x)
  
    xcircularp <- attr(as.circular(xx[,1]), "circularp")
    type <- xcircularp$type
    units <- xcircularp$units
    template <- xcircularp$template
    modulo <- xcircularp$modulo
    zero <- xcircularp$zero
    rotation <- xcircularp$rotation

    xlim <- shrink * xlim
    ylim <- shrink * ylim
    midx <- 0.5 * (xlim[2] + xlim[1])
    xlim <- midx + (1 + tol) * 0.5 * c(-1, 1) * (xlim[2] - xlim[1])
    midy <- 0.5 * (ylim[2] + ylim[1])
    ylim <- midy + (1 + tol) * 0.5 * c(-1, 1) * (ylim[2] - ylim[1])
    oldpin <- par("pin")
    xuin <- oxuin <- oldpin[1]/diff(xlim)
    yuin <- oyuin <- oldpin[2]/diff(ylim)
    if (missing(uin)) {
        if (yuin > xuin) yuin <- xuin
        else xuin <- yuin
    } else {
        if (length(uin) == 1) uin <- uin * c(1, 1)
        if (any(c(xuin, yuin) < uin)) stop("uin is too large to fit plot in")
        xuin <- uin[1]; yuin <- uin[2]
    }    
    xlim <- midx + oxuin/xuin * c(-1, 1) * diff(xlim) * 0.5
    ylim <- midy + oyuin/yuin * c(-1, 1) * diff(ylim) * 0.5
    plot(cos(seq(0, 2 * pi, length = 1000)), sin(seq(0, 2 * pi, length = 1000)), axes = FALSE, xlab = "", ylab = "", main = main, type = "l", xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
    
    if (missing(bins)) {
         bins <- NROW(x)
    } else {
         bins <- round(bins)
         if (bins<=0) stop("bins must be non negative")
    }

    if (axes) {
        axis.circular(units = units, template=template, zero=zero, rotation=rotation, digits=digits, cex=cex, tcl=tcl, tcl.text=tcl.text)
    }

    if (missing(col)) {
    col <- seq(nseries)
    } else {
    if (length(col)!=nseries) {
        col <- rep(col, nseries)[1:nseries]
    }
    }
    pch <- rep(pch, nseries, length.out=nseries)
    
    if (!is.logical(ticks)) stop("ticks must be logical")

    arc <- (2 * pi)/bins
    pos.bins <- ((1:nseries)-1/2)*arc/nseries-arc/2
            
    if (ticks) {
        at <- (0:bins)/bins*2*pi
        if (rotation=="clock") at <- -at
        at <- at + zero

        ticks.circular(circular(x=at, type="angles", units="radians", modulo="asis", zero=zero, rotation=rotation), tcl=tcl)
    }

    for (iseries in 1:nseries) {
      
     x <- xx[,iseries]
         x <- as.circular(x)    
         x <- conversion.circular(x, units="radians")
     n <- length(x)

     if (!stack) {
             if (rotation=="clock") x <- -x
             x <- x + zero 
             z <- cos(x)
             y <- sin(x)
         r <- 1+(iseries-1)*sep*shrink
         points.default(z*r, y*r, cex=cex, pch=pch[iseries], col = col[iseries], ...)
     } else {
         bins.count <- c(1:bins)
         xmod <- x %% (2 * pi)
         for (i in 1:bins) {
          bins.count[i] <- sum(xmod <= i * arc & xmod > (i - 1) * arc)
         }
         mids <- seq(arc/2, 2 * pi - pi/bins, length = bins) + pos.bins[iseries]
             if (rotation=="clock") mids <- -mids
             mids <- mids + zero
 
         index <- cex*sep
         for (i in 1:bins) {
          if (bins.count[i] != 0) {
              for (j in 0:(bins.count[i] - 1)) {
               r <- 1 + j * index
               z <- r * cos(mids[i])
               y <- r * sin(mids[i])
               points.default(z, y, cex=cex, pch=pch[iseries], col=col[iseries], ...)
              }
           }
          }
     }
    }
return(invisible(list(zero=zero, rotation=rotation, next.points=nseries*sep)))
}

#############################################################
#                                                           #
#   points.circular function                                #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: September, 22, 2003                               #
#   Version: 0.1-2                                          #
#                                                           #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################
 
points.circular <- function(x, pch = 16, cex = 1, stack = FALSE, sep = 0.025, shrink=1, bins, col, next.points, plot.info, zero, rotation, ...) {

    if (is.matrix(x) | is.data.frame(x)) {
        nseries <- ncol(x)
    } else {
        nseries <- 1
    }
    xx <- as.data.frame(x)
  
    xcircularp <- attr(as.circular(xx[,1]), "circularp")
    type <- xcircularp$type
    modulo <- xcircularp$modulo
    if (missing(plot.info)) {
        if (missing(zero)) zero <- xcircularp$zero
        if (missing(rotation)) rotation <- xcircularp$rotation
        if (missing(next.points)) next.points <- 0
    } else {
        zero <- plot.info$zero
        rotation <- plot.info$rotation
        if (missing(next.points))
            next.points <- plot.info$next.points
    }
        
    x <- conversion.circular(x, units="radians")

    if (missing(bins)) {
    bins <- NROW(x)
    } else {
    bins <- round(bins)
    if (bins<=0) stop("bins must be non negative")
    }

    if (is.matrix(x) || is.data.frame(x)) {
        nseries <- ncol(x)
    } else {
        nseries <- 1
    }
 
    if (missing(col)) {
    col <- seq(nseries)
    } else {
    if (length(col)!=nseries) {
        col <- rep(col, nseries)[1:nseries]
    }
    }
    pch <- rep(pch, nseries, length.out=nseries)

    arc <- (2 * pi)/bins
    pos.bins <- ((1:nseries)-1/2)*arc/nseries-arc/2
            
    for (iseries in 1:nseries) {
     x <- xx[,iseries] 
         x <- as.circular(x)
         x <- conversion.circular(x, units="radians")
         n <- length(x)

     if (!stack) {
             if (rotation=="clock") x <- -x
             x <- x + zero
             z <- cos(x)
             y <- sin(x)
         r <- 1+next.points+(iseries-1)*sep*shrink
         points.default(z*r, y*r, cex=cex, pch=pch[iseries], col = col[iseries], ...)
     } else {
         bins.count <- c(1:bins)

         for (i in 1:bins) {
          bins.count[i] <- sum(x <= i * arc & x > (i - 1) * arc)
         }
         mids <- seq(arc/2, 2 * pi - pi/bins, length = bins) + pos.bins[iseries]
             if (rotation=="clock") mids <- -mids
             mids <- mids + zero
 
         index <- cex*sep
         for (i in 1:bins) {
          if (bins.count[i] != 0) {
              for (j in 0:(bins.count[i] - 1)) {
               r <- 1 + j * index
               z <- r * cos(mids[i])
               y <- r * sin(mids[i])
               points.default(z, y, cex=cex, pch=pch[iseries], col=col[iseries], ...)
              }
           }
          }
     }
    }
return(invisible(list(zero=zero, rotation=rotation, next.points=next.points+nseries*sep)))
}
