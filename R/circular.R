#############################################################
#                                                           #
#   circular function                                       #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: September, 22, 2003                               #
#   Version: 0.6-1                                          #
#                                                           #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################

circular <- function(x, type=c("angles", "directions"), units=c("radians", "degrees"), template=c("none", "geographics"), modulo=c("asis", "2pi", "pi"), zero=0, rotation=c("counter", "clock"), names) {

    type <- match.arg(type)
    units <- match.arg(units)
    template <- match.arg(template) 
    modulo <- match.arg(modulo)
    rotation <- match.arg(rotation)
  
    if (template=="geographics") {
        zero <- pi/2
        rotation <- "clock"
    }

    if (is.data.frame(x)) x <- as.matrix(x)

    if (is.matrix(x)) {
    nseries <- ncol(x)
    ndata <- nrow(x)
    if (missing(names)) {
            names <- if(!is.null(dimnames(x))) colnames(x) else paste("Circular", seq(nseries), sep="")
    }
        dimnames(x) <- list(NULL, names)
    } else {
        nseries <- 1
    ndata <- length(x)
    }
####    if (ndata == 0) stop("circular object must have one or more observations")
    if (modulo!="asis") {
    if (modulo=="2pi") {
        ang <- 2
    } else {
        ang <- 1
    }
    if (units=="radians") {
        x <- x %% (ang*pi)
    } else {
        x <- x %% (ang*180)
    }
    }

    attr(x, "circularp") <- list(type=type, units=units, template=template, modulo=modulo, zero=zero, rotation=rotation) #-- order is fixed
    attr(x, "class") <- "circular"
    return(x)
}

#############################################################
#                                                           #
#       conversion.circular function                        #
#   Author: Claudio Agostinelli                         #
#   E-mail: claudio@unive.it                            #
#   Date: July, 31, 2003                                #
#   Version: 0.1-1                                      #
#                                                           #
#   Copyright (C) 2003 Claudio Agostinelli              #
#                                                           #
#############################################################

conversion.circular <- function(x, units=c("radians", "degrees")) {
    units <- match.arg(units)
    x <- as.circular(x)
    value <- attr(x, "circularp")
    unitsp <- value$units
          
    if (unitsp=="degrees" & units=="radians") {
    x <- x/180*pi
    } else if (unitsp=="radians" & units=="degrees") {
               x <- x/pi*180
    }
    value$units <- units
    circularp(x) <- value 
    
    return(x)
}

#############################################################
#                                                           #
#   circularp function                                  #
#   Author: Claudio Agostinelli                         #
#   E-mail: claudio@unive.it                            #
#   Date: March, 7, 2003                                #
#   Version: 0.1                                        #
#                                                           #
#   Copyright (C) 2003 Claudio Agostinelli              #
#                                                           #
#############################################################
 
circularp <- function(x) attr(x, "circularp")

#############################################################
#                                                           #
#   circularp<- function                                    #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: November, 18, 2003                                #
#   Version: 0.2-2                                          #
#                                                           #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################
 
"circularp<-" <- function(x, value) {
    cl <- class(x)
    if (length(value)!=6) stop("value must have six elements")

    type <- value$type
    units <- value$units
    template <- value$template
    modulo <- value$modulo
    zero <- value$zero
    rotation <- value$rotation
 
    if (type!="angles" & type!="directions") stop("type (value[1]) must be 'angles', 'directions' or 'geographics'")

    if (units!="radians" & units!="degrees") stop("units (value[2]) must be 'radians' or 'degrees'")

    if (template!="none" & template!="geographics") stop("template (value[3]) must be 'none' or 'geographics'")
    
    if (modulo!="asis" & modulo!="2pi" &  modulo!="pi") stop("modulo (value[4]) must be 'asis' or 'pi' or '2pi'")

    if (rotation!="clock" & rotation!="counter") stop("rotation (value[6]) must be 'clock' or 'counter'")

    attr(x, "circularp") <- value
    if (inherits(x, "circular") && is.null(value))
        class(x) <- cl["circular" != cl]
    return(x)
}

#############################################################
#                                                           #
#   is.circular function                                #
#   Author: Claudio Agostinelli                         #
#   E-mail: claudio@unive.it                            #
#   Date: March, 7, 2003                                #
#   Version: 0.1                                        #
#                                                           #
#   Copyright (C) 2003 Claudio Agostinelli              #
#                                                           #
#############################################################
 
is.circular <- function (x) inherits(x, "circular")

#############################################################
#                                                           #
#   [.circular function                                 #
#   Author: Claudio Agostinelli                         #
#   E-mail: claudio@unive.it                            #
#   Date: March, 7, 2003                                #
#   Version: 0.1                                        #
#                                                           #
#   Copyright (C) 2003 Claudio Agostinelli              #
#                                                           #
#############################################################
 
"[.circular" <- function(x, i, ...) {
    y <- NextMethod("[", ...)
    class(y) <- class(x)
    attr(y, "circularp") <- attr(x, "circularp")
    return(y)
}

#############################################################
#                                                           #
#   print.circular function                             #
#   Author: Claudio Agostinelli                         #
#   E-mail: claudio@unive.it                            #
#   Date: June, 21, 2003                                #
#   Version: 0.2                                        #
#                                                           #
#   Copyright (C) 2003 Claudio Agostinelli              #
#                                                           #
#############################################################
 
print.circular <- function(x, info=TRUE, ...) {
    x.orig <- x
    x <- as.circular(x)
    if (info) {
    xcircularp <- attr(x, "circularp")
    type <- xcircularp$type
    units <- xcircularp$units
        template <- xcircularp$template
    modulo <- xcircularp$modulo
        zero <- xcircularp$zero
        rotation <- xcircularp$rotation

        cat("Circular Data: \nType =", type,
               "\nUnits =", units,
               "\nTemplate =", template, 
               "\nModulo =", modulo,
               "\nZero =", zero,
               "\nRotation =", rotation, "\n")
    }       
    attr(x, "class") <- attr(x, "circularp") <- attr(x, "na.action") <- NULL
    NextMethod("print", x, quote = FALSE, right = TRUE, ...)
    invisible(x.orig)
}
