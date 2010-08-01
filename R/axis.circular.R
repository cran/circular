#############################################################
#                                                           #
#   axis.circular function                                  #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: October, 19, 2009                                 #
#   Version: 0.4-3                                          #
#                                                           #
#   Copyright (C) 2009 Claudio Agostinelli                  #
#                                                           #
#############################################################
 
axis.circular <- function(at=NULL, labels=NULL,  units = NULL, template=NULL, modulo=NULL, zero=NULL, rotation=NULL, tick=TRUE, lty, lwd, cex, col, font, tcl=0.025, tcl.text=0.125, digits=2) {

   if (missing(cex)) cex <- par("cex.axis")
   if (missing(col)) col <- par("col.axis")
   if (missing(font)) font <- par("font.axis")
   if (missing(lty)) lty <- par("lty")
   if (missing(lwd)) lwd <- par("lwd")

   if (is.null(at)) {
      if (is.null(template) | template=="none" | template=="geographics") {
        at <- circular(c(0, pi/2, pi, 3/2*pi))
      } else if (template=="clock24") {
        at <- circular(seq(0, 23), units="hours")
      } else if (template=="clock12") {
        at <- circular(seq(0, 11), units="hours", modulo="pi")
      }
   }
   at <- na.omit(at)
   atcircularp <- attr(as.circular(at), "circularp")
   type <- atcircularp$type
   if (is.null(modulo))
      modulo <- atcircularp$modulo
   if (is.null(units)) 
      units <- atcircularp$units
   if (is.null(template))
      template <- atcircularp$template
   if (template=="geographics" | template=="clock24") {
      zero <- pi/2
      rotation <- "clock"
   } else if (template=="clock12") {
      zero <- pi/2
      rotation <- "clock"
      modulo <- "pi"
   } else {
      if (is.null(zero))
         zero <- atcircularp$zero
      if (is.null(rotation))
         rotation <- atcircularp$rotation
   }

   atasis <- at
   attr(atasis, "circularp") <- attr(atasis, "class") <- NULL
   attext <- atasis/pi
   if (modulo=="2pi") {
      if (units=="radians") {
         atasis <- atasis%%(2*pi)
      } else if (units=="degrees") {
         atasis <- atasis%%(360)
      } else if (units=="hours") {
         atasis <- atasis%%(24)
      }
      attext <- attext%%2
    } else if (modulo=="pi") {
              if (units=="radians") {
                 atasis <- atasis%%pi
              } else if (units=="degrees") {
                 atasis <- atasis%%180
              } else if (units=="hours") {
                 atasis <- atasis%%12
              }
              attext <- attext%%1
          }
   attext <- round(attext, digits=digits)
   
   if (template=="clock12")
     at <- 2*at
   at <- conversion.circular(at, units="radians", modulo="2pi", zero=0, rotation="counter")
   attr(at, "circularp") <- attr(at, "class") <- NULL
   if (rotation=="clock")
            at <- -at
   at <- at+zero
   if (is.null(labels)) {
      if (length(atasis)==4 && all(atasis==c(0, pi/2, pi, 3/2*pi))) { 
         if (template=="geographics") {
            labels <- c("N", "E", "S", "W")         
         } else {
            if (units=="radians") {
               labels <- c("0", expression(frac(pi,2)), expression(pi), expression(frac(3*pi,2)))
            } else if (units=="degrees") {
               labels <- c("0", "90", "180", "270")      
            }
         }
      } else if (units=="degrees") {
                labels <- as.character(round(atasis, digits=digits))
      } else if (units=="hours") {
                labels <- as.character(round(atasis, digits=digits))
      }
  }

  if (!is.null(labels) && length(at)!=length(labels))
     stop("'at' and 'labels' must have the same length")
  
   AxisCircularRad(at, units, labels, attext, tick, tcl, tcl.text, cex, col, lty, lwd)
}

AxisCircularRad <- function(at, units, labels, attext, tick, tcl, tcl.text, cex, col, lty, lwd) {
#### at must be in radians, counter, zero=0
   r <- 1+tcl*c(-1/2,1/2)
   r.l <- 1-tcl.text 
   z <- cos(at)
   y <- sin(at)
  
   for (i in 1:length(at)) {
      if (tick) {
         lines.default(z[i]*r, y[i]*r, col=col, lty=lty, lwd=lwd)
      }
      if (is.null(labels)) {
          labeltext <- substitute(at*pi, list(at=attext[i]))
      } else {
          labeltext <- labels[i]
      }
      text.default(z[i]*r.l, y[i]*r.l, labeltext, cex=cex, col=col)    
   }  
   text(0, 0, "+", cex=1)
}

