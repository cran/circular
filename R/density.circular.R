#############################################################
#                                                           #
#   density.circular function                               #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: November, 19, 2003                                #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1-5                                           #
#                                                           #
#############################################################

density.circular <- function(x, z, bw, adjust = 1, type = c("K", "L"), kernel = c("vonmises", "wrappednormal"), na.rm = FALSE, from=0, to=2*pi, n=512, K=10, ...) {

    name <- deparse(substitute(x))
    xx <- x
    x <- as.circular(x)

    xcircularp <- circularp(x)
    type <- xcircularp$type
    units <- xcircularp$units
    template <- xcircularp$template
    zero <- xcircularp$zero
    rotation <- xcircularp$rotation
    
    x <- conversion.circular(x, units="radians")
  
    kernel <- match.arg(kernel)

    if (!is.numeric(n))
        stop("argument must be numeric")
    n <- round(n)
    if (n <=0)
         stop("argument must be integer and positive")     

    if (!is.numeric(from))
        stop("argument must be numeric")      
    if (!is.numeric(to))
        stop("argument must be numeric")      
    if (!is.finite(from)) 
        stop("non-finite `from'")
    if (!is.finite(to)) 
        stop("non-finite `to'")
    
    if (!is.numeric(x)) 
        stop("argument must be numeric")
    x <- as.vector(x)
    x.na <- is.na(x)
    if (any(x.na)) {
        if (na.rm) 
            x <- x[!x.na]
        else stop("x contains missing values")
    }
    
    nx <- length(x)
    x.finite <- is.finite(x)
    if (any(!x.finite)) {
        x <- x[x.finite]
        nx <- sum(x.finite)
    }

    if (missing(z)) {
        z <- zz <- circular(seq(from=from, to=to, length=n), units="radians", type=type, template=template, zero=zero, rotation=rotation)
        
    } else {
        z <- zz <- as.circular(z)
        z <- conversion.circular(z, units="radians")
        if (!is.numeric(z))
            stop("argument 'z' must be numeric")
        namez <- deparse(substitute(z))
        z <- as.vector(z)
        z.na <- is.na(z)
        if (any(z.na)) {
            if (na.rm) {
                z <- z[!z.na]
            } else {
                stop("z contains missing values")
            }
        }
    
        nz <- length(z)
        z.finite <- is.finite(z)
        if (any(!z.finite)) {
            z <- z[z.finite]
            nz <- sum(z.finite)
        }
    }

    bw <- adjust * bw
    if (!is.numeric(bw))
        stop("argument must be numeric")        
    if (!is.finite(bw)) 
        stop("non-finite `bw'")
    if (bw <= 0) 
        stop("`bw' is not positive.")

    if (kernel=="vonmises") {
        y <- sapply(z, dvonmises, mu=x, kappa=bw)
    } else if (kernel=="wrappednormal") {
        y <- sapply(z, dwrappednormal, mu=x, sd=bw, K=K)
    } else {
        stop("other kernels not implemented yet")
    }
    y <- apply(y, 2, sum)/nx

    if (units=="degrees") xx <- conversion.circular(xx, units="degrees")

    structure(list(data = xx, x = zz, y = y, bw = bw, n = nx, kernel=kernel, call = match.call(), data.name=name, has.na = FALSE), class = "density.circular")
} 

#############################################################
#                                                           #
#   plot.density.circular function                          #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: August, 01, 2003                                  #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#   Version 0.2-1                                           #
#                                                           #
#############################################################

plot.density.circular <- function(x, main = NULL, xlab = NULL, ylab = "Density circular", type = "l", zero.line = TRUE, points.plot=FALSE, points.col=1, points.pch=1, plot.type = c("circle", "line"), axes=TRUE, ticks=TRUE, bins, shrink=1, tcl=0.025, tol = 0.04, uin, xlim=c(-1, 1), ylim=c(-1, 1), ...) {

    x$x <- conversion.circular(x$x, units="radians")
    x$data <- conversion.circular(x$data, units="radians")
  
    plot.type <- match.arg(plot.type)
    if (missing(bins)) {
	bins <- NROW(x)
    } else {
	bins <- round(bins)
	if (bins<=0) stop("bins must be non negative")
    }
    
    if (is.null(xlab)) 
        xlab <- paste("N =", x$n, "  Bandwidth =", formatC(x$bw))
    if (is.null(main)) 
        main <- deparse(x$call)

    if (plot.type == "line") {
        xorder <- order(x$x)
        x$x <- x$x[xorder]
        x$y <- x$y[xorder] 
        plot.default(x, main = main, xlab = xlab, ylab = ylab, type = type, ...)
        if (zero.line) 
            abline(h = 0, lwd = 0.1, col = "gray")
        if (points.plot)
            points(x$data, rep(min(x$y),length(x$data)), col=points.col, pch=points.pch)
    } else {
        x$x <- as.circular(x$x)
        xcircularp <- attr(x$x, "circularp")
        xtype <- xcircularp$type
        units <- xcircularp$units
        template <- xcircularp$template
        modulo <- xcircularp$modulo
        zero <- xcircularp$zero
        rotation <- xcircularp$rotation
    
        x$x <- conversion.circular(x$x, units="radians")

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
       plot(cos(seq(0, 2 * pi, length = 1000)), sin(seq(0, 2 * pi, length = 1000)), axes = FALSE, xlab = "", ylab = "", main = "", type = "l", xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")
 
        if (rotation=="clock") x$x <- -x$x
        x$x <- x$x + zero
            
        if (axes) {
	    axis.circular(units = units, template=template, modulo = modulo, zero=zero, rotation=rotation)
        }
 
        if (ticks) {
            at <- (0:bins)/bins*2*pi
            if (rotation=="clock") at <- -at
            at <- at + zero

            ticks.circular(circular(x=at, type="angles", units="radians", modulo="asis", zero=zero, rotation=rotation), tcl=tcl)
        }

        z <- (x$y+1)*cos(x$x)
        y <- (x$y+1)*sin(x$x)
        xorder <- order(x$x)
        z <- z[xorder]
        y <- y[xorder] 
        lines(x=z, y=y, type = type, ...)
        if (points.plot) {
            points.circular(x$data, col=points.col, pch=points.pch)
        }
    }
}

#############################################################
#                                                           #
#   lines.density.circular function                         #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: July, 23, 2003                                    #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1                                             #
#                                                           #
#############################################################

lines.density.circular <- function(x, type = "l", zero.line = TRUE, points.plot=FALSE, points.col=1, points.pch=1, plot.type = c("circle", "line"), bins, shrink=1, tcl=0.025, ...) {

    x$x <- conversion.circular(x$x, units="radians")
    x$data <- conversion.circular(x$data, units="radians")
  
    plot.type <- match.arg(plot.type)
    if (missing(bins)) {
	bins <- NROW(x)
    } else {
	bins <- round(bins)
	if (bins<=0) stop("bins must be non negative")
    }
    
    if (plot.type == "line") {
        xorder <- order(x$x)
        x$x <- x$x[xorder]
        x$y <- x$y[xorder] 
        lines.default(x, type = type, ...)
        if (zero.line) 
            abline(h = 0, lwd = 0.1, col = "gray")
        if (points.plot)
            points(x$data, rep(min(x$y),length(x$data)), col=points.col, pch=points.pch)
    } else {
        x$x <- as.circular(x$x)
        xcircularp <- attr(x$x, "circularp")
        xtype <- xcircularp$type
        units <- xcircularp$units
        template <- xcircularp$template
        modulo <- xcircularp$modulo
        zero <- xcircularp$zero
        rotation <- xcircularp$rotation
    
        x$x <- conversion.circular(x$x, units="radians")

        if (rotation=="clock") x$x <- -x$x
        x$x <- x$x + zero
            
        z <- (x$y+1)*cos(x$x)
        y <- (x$y+1)*sin(x$x)
        xorder <- order(x$x)
        z <- z[xorder]
        y <- y[xorder] 
        lines(x=z, y=y, type = type, ...)
        if (points.plot) {
            points.circular(x$data, col=points.col, pch=points.pch)
        }
    }
}


#############################################################
#                                                           #
#   print.density.circular function                         #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: July, 23, 2003                                    #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1                                             #
#                                                           #
#############################################################

print.density.circular <- function(x, digits=NULL, ...)
{
    cat("\nCall:\n\t",deparse(x$call),
	"\n\nData: ",x$data.name," (",x$n," obs.);",
	"\tBandwidth 'bw' = ",formatC(x$bw,digits=digits), "\n\n",sep="")
    print(summary(as.data.frame(x[c("x","y")])), digits=digits, ...)
    invisible(x)
}

