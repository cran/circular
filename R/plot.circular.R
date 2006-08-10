#############################################################
#                                                           #
#   plot.circular function                                  #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: June, 07, 2006                                    #
#   Version: 0.3-1                                          #
#                                                           #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#############################################################
 
plot.circular <- function(x, pch=16, cex=1, stack=FALSE, axes=TRUE, sep=0.025, shrink=1, bins=NULL, ticks=FALSE, tcl=0.025, tcl.text=0.125, col=NULL, tol=0.04, uin=NULL, xlim=c(-1, 1), ylim=c(-1, 1), digits=2, units=NULL, template=NULL, zero=NULL, rotation=NULL, main="", xlab="", ylab="", ...) {

   if (is.matrix(x) | is.data.frame(x)) {
      nseries <- ncol(x)
   } else {
      nseries <- 1
   }
   xx <- as.data.frame(x)
  
   xcircularp <- attr(as.circular(xx[,1]), "circularp")
   type <- xcircularp$type
   modulo <- xcircularp$modulo
   if (is.null(units)) 
      units <- xcircularp$units
   if (is.null(template))
      template <- xcircularp$template
   if (template=="geographics") {
      zero <- pi/2
      rotation <- "clock"
   } else {
      if (is.null(zero))
         zero <- xcircularp$zero
      if (is.null(rotation))
         rotation <- xcircularp$rotation
   }
   
   CirclePlotRad(xlim, ylim, uin, shrink, tol, 1000, main=main, xlab=xlab, ylab=ylab)
    
   if (is.null(bins)) {
      bins <- NROW(x)
   } else {
      bins <- round(bins)
      if (bins<=0)
         stop("bins must be non negative")
   }

   if (is.null(col)) {
      col <- seq(nseries)
   } else {
      if (length(col)!=nseries) {
         col <- rep(col, nseries)[1:nseries]
      }
   }
   pch <- rep(pch, nseries, length.out=nseries)

   if (axes) {
      axis.circular(units=units, template=template, zero=zero, rotation=rotation, digits=digits, cex=cex, tcl=tcl, tcl.text=tcl.text)
   }
   
   if (!is.logical(ticks))
      stop("ticks must be logical")
   
   if (ticks) {
      at <- circular((0:bins)/bins*2*pi, zero=zero, rotation=rotation)
      ticks.circular(at, tcl=tcl)
   }

   for (iseries in 1:nseries) {      
      x <- xx[,iseries]
      x <- na.omit(x)
      n <- length(x)      
      if (n) {
         x <- conversion.circular(x, units="radians")
         attr(x, "circularp") <- attr(x, "class") <- NULL
         if (rotation=="clock")
            x <- -x
         x <- x+zero
         x <- x%%(2*pi)
         PointsCircularRad(x, bins, stack, col, pch, iseries, nseries, sep, 0, shrink, cex, ...) 
      }
   }
return(invisible(list(zero=zero, rotation=rotation, next.points=nseries*sep)))
}

CirclePlotRad <- function(xlim=c(-1,1), ylim=c(-1,1), uin=NULL, shrink=1, tol=0.04, n=1000, ...) {
   xlim <- shrink * xlim
   ylim <- shrink * ylim
   midx <- 0.5 * (xlim[2] + xlim[1])
   xlim <- midx + (1 + tol) * 0.5 * c(-1, 1) * (xlim[2] - xlim[1])
   midy <- 0.5 * (ylim[2] + ylim[1])
   ylim <- midy + (1 + tol) * 0.5 * c(-1, 1) * (ylim[2] - ylim[1])
   oldpin <- par("pin")
   xuin <- oxuin <- oldpin[1]/diff(xlim)
   yuin <- oyuin <- oldpin[2]/diff(ylim)
   if (is.null(uin)) {
      if (yuin > xuin)
         yuin <- xuin
      else
         xuin <- yuin
   } else {
      if (length(uin) == 1)
         uin <- uin * c(1, 1)
      if (any(c(xuin, yuin) < uin))
         stop("uin is too large to fit plot in")
      xuin <- uin[1]; yuin <- uin[2]
   }
   xlim <- midx + oxuin/xuin * c(-1, 1) * diff(xlim) * 0.5
   ylim <- midy + oyuin/yuin * c(-1, 1) * diff(ylim) * 0.5
   plot.default(cos(seq(0, 2 * pi, length = n)), sin(seq(0, 2 * pi, length = n)), axes = FALSE, type = "l", xlim=xlim, ylim=ylim, xaxs="i", yaxs="i", ...)
}

