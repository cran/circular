
###############################################################
#                                                             #
#       Original Splus: Ulric Lund                            #
#       E-mail: ulund@calpoly.edu                             #
#                                                             #
###############################################################

#############################################################
#                                                           #
#   rose.diag function                                      #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: April, 11, 2005                                   #
#   Copyright (C) 2005 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1-4                                           #
#                                                           #
#############################################################

rose.diag <- function(x, pch = 16, axes = TRUE, shrink = 1, bins, ticks = TRUE, tcl=0.025, col, tol = 0.04, uin, xlim=c(-1, 1), ylim=c(-1, 1), prop = 1, main=NULL, ...) {
  
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
    axis.circular(units = units, template=template, modulo = modulo, zero=zero, rotation=rotation)
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
# Add to remove NA values from each series
        x <- na.omit(x)
        x <- as.circular(x)    
        x <- conversion.circular(x, units="radians")
        if (rotation=="clock") x <- -x
        x <- x + zero 
        x <- x %% (2 * pi)
        n <- length(x)
    freq <- c(1:bins)
    arc <- (2 * pi)/bins
    for(i in 1:bins) {
        freq[i] <- sum(x <= i * arc & x > (i - 1) * arc)
    }
    rel.freq <- freq/n
    radius <- sqrt(rel.freq) * prop
    sector <- seq(0, 2 * pi - (2 * pi)/bins, length = bins)
    mids <- seq(arc/2, 2 * pi - pi/bins, length = bins)
    for(i in 1:bins) {
        if(rel.freq[i] != 0) {
            lines.default(c(0, radius[i] * cos(sector[i])), c(0, radius[i] * sin(sector[i])), col=col[iseries], ...)
            lines.default(c(0, radius[i] * cos(sector[i] + (2 * pi)/bins)), c(0, radius[i] * sin(sector[i] + (2 * pi)/bins)), col=col[iseries], ...)
            lines.default(c(radius[i] * cos(sector[i]), radius[i] * cos(sector[i] + (2 * pi)/bins)), c(radius[i] * sin(sector[i]), radius[i] * sin(sector[i] + (2 * pi)/bins)), col=col[iseries], ...)
        }
    }
   }
return(invisible(list(zero=zero, rotation=rotation, next.points=0)))    
}
