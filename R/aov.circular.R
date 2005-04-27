
###############################################################
#                                                             #
#       Original Splus: Ulric Lund                            #
#       E-mail: ulund@calpoly.edu                             #
#                                                             #
###############################################################

#############################################################
#                                                           #
#   aov.circular function                                   #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: April, 13, 2005                                   #
#   Version: 0.1-2                                          #
#                                                           #
#   Copyright (C) 2005 Claudio Agostinelli                  #
#                                                           #
#############################################################

aov.circular <- function(x, group, kappa=NULL, method=c("F.test", "LRT"), F.mod=TRUE){

    # Handling missing values
    ok <- complete.cases(x, group)
    x <- x[ok]
    group <- group[ok]
    if (length(x)==0 | length(table(group)) < 2) {
        warning("No observations or no groups (at least after removing missing values)")
        return(NULL)
    }
    
    x <- as.circular(x)
    xcircularp <- circularp(x)
    units <- xcircularp$units
    x <- conversion.circular(x, units="radians")
    x <- x%%(2*pi)
    method <- match.arg(method)
    ns        <- tapply(x, group, FUN=length)
    resultant <- tapply(x, group, FUN=function(x) rho.circular(x)*length(x))
    mean.dirs <- tapply(x, group, FUN=mean.circular)
    kappas    <- tapply(x, group, FUN=function(x) mle.vonmises(x)$kappa)
    grps <- length(resultant)
    n <- length(group)
    res.all <- rho.circular(x)*n
    mean.dir.all <- mean.circular(x)
    kappa.all <- mle.vonmises(x)$kappa

    if (method=="F.test"){
        if (!is.null(kappa))
            warning("Specified value of kappa is not used in the F-test")
        sum.res <- sum(resultant)
        df <- c(grps-1, n-grps, n-1)
        SS <- c(sum.res - res.all, n-sum.res, n-res.all) 
        MS <- SS/df
        if (F.mod==TRUE) {
            F.stat <- (1+3/(8*kappa.all))*MS[1]/MS[2]
        } else {
            F.stat <- MS[1]/MS[2]
        }
        p.value <- 1-pf(F.stat, grps-1,n-grps)
    } else {
        if (is.null(kappa))
            kappa <- kappa.all
        stat1 <- 1-1/(4*kappa)*A1(kappa)*(sum(1/ns)-1/n)
        stat2 <- 2*kappa*sum(resultant*(1-cos(mean.dirs-mean.dir.all)))
        chisq.stat <- stat1*stat2
        p.value <- 1-pchisq(chisq.stat, grps-1)
       
    }
    mean.dir.all <- conversion.circular(mean.dir.all, units=units)
    if (units=="degrees") {
        mean.dirs <- mean.dirs/pi*180
        
    }
    
    attr(mean.dirs, "circularp") <- xcircularp
    attr(mean.dirs, "class") <- "circular"
 
    result <- list()
    result$call <- match.call()
    result$mu <- mean.dirs 
    result$mu.all <- mean.dir.all
    result$kappa <- kappas 
    result$kappa.all <- kappa.all
    result$rho <- resultant
    result$rho.all <- res.all
    result$method <- method

    if (method=="F.test") { 
        result$df <- df 
        result$SSE <- SS
        result$MSE <- MS
        result$statistic <- F.stat
    } else {
        result$df <- grps-1
        result$statistic <- chisq.stat
    }

    result$p.value <- p.value
    class(result) <- "aov.circular"
    return(result)
}

#############################################################
#                                                           #
#   print.aov.circular function                             #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: April, 13, 2005                                   #
#   Version: 0.2                                            #
#                                                           #
#   Copyright (C) 2005 Claudio Agostinelli                  #
#                                                           #
#############################################################

print.aov.circular <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="") 
 
    if (x$method=="F.test") {
        result.matrix <- cbind(x$df, x$SSE, x$MSE, c(x$statistic,NA,NA), c(x$p.value,NA,NA))
        dimnames(result.matrix) <- list(c("Between","Within","Total"),c("df", "SSE", "MSE", "F", "p"))
        cat("\n", "Circular Analysis of Variance: High Concentration F-Test", "\n", "\n")
        print(result.matrix, digits=digits)
        cat("\n \n")
    } else {
        cat("\n", "Circular Analysis of Variance: Likelihood Ratio Test", "\n", "\n")
        cat(" df:     ", format(x$df, digits=digits), "\n ChiSq:  ", format(x$statistic, digits=digits), "\n p.value:",  format(x$p.value, digits=digits), "\n \n")
    }
    invisible(x)
}

