
###############################################################
#                                                             #
#       Original Splus: Ulric Lund                            #
#       E-mail: ulund@calpoly.edu                             #
#                                                             #
###############################################################

#############################################################
#                                                           #
#   rao.test function                                       #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: July, 25, 2003                                    #
#   Version: 0.1                                            #
#                                                           #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################

rao.test <- function(..., alpha = 0) {
        y <- list(...)
        x <- list()
        for (i in 1:length(y)) {
             if (is.data.frame(y[[i]])) {
                 x <- c(x, as.list(y[[i]]))
             } else if (is.matrix(y[[i]])) {
                        for (j in 1:ncol(y[[i]])) {
                             x <- c(x, list(y[[i]][,j]))
                        }
             } else if (is.list(y[[i]])) {
                        x <- c(x, y[[i]])
             } else {
                 x <- c(x, list(y[[i]]))
             }
        }
        if (length(x)<2) stop("There must be at least two samples")
        for (i in 1:length(x)) {
             x[[i]] <- as.circular(x[[i]])
             x[[i]] <- conversion.circular(x[[i]], units="radians")
             attr(x[[i]], "circularp") <- attr(x[[i]], "class") <- NULL
        }
        if (!any(c(0, 0.01, 0.025, 0.05, 0.1, 0.15)==alpha)) stop("'alpha' must be one of the following values: 0, 0.01, 0.025, 0.05, 0.1, 0.15")
    n <- unlist(lapply(x, length))
    k <- length(x)
    c.data <- lapply(x, cos)
    s.data <- lapply(x, sin)
    x <- unlist(lapply(c.data, mean))
    y <- unlist(lapply(s.data, mean))
    s.co <- unlist(lapply(c.data, var))
    s.ss <- unlist(lapply(s.data, var))
    s.cs <- c(1:k)
    for(i in 1:k) {
        s.cs[i] <- var(c.data[[i]], s.data[[i]])
    }
    s.polar <- 1/n * (s.ss/x^2 + (y^2 * s.co)/x^4 - (2 * y * s.cs)/x^3)
    tan <- y/x
    H.polar <- sum(tan^2/s.polar) - (sum(tan/s.polar))^2/sum(1/s.polar)
    U <- x^2 + y^2
    s.disp <- 4/n * (x^2 * s.co + y^2 * s.ss + 2 * x * y * s.cs)
    H.disp <- sum(U^2/s.disp) - (sum(U/s.disp))^2/sum(1/s.disp)
        
    result <- list(statistic=c(H.polar, H.disp), df=k-1, p.value=c((1 - pchisq(H.polar, k - 1)), (1 - pchisq(H.disp, k - 1))), alpha=alpha)
    class(result) <- "rao.test"
    return(result)
}

#############################################################
#                                                           #
#   print.rao.test function                                 #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: November, 19, 2003                                #
#   Version: 0.1-1                                          #
#                                                           #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################

print.rao.test <- function(x, digits=4, ...) {
        statistic <- x$statistic
        p.value <- x$p.value
        alpha <- x$alpha
        df <- x$df
    cat("\n")
    cat("Rao's Tests for Homogeneity", "\n")
    if(alpha == 0) {
        cat("\n")
        cat("       Test for Equality of Polar Vectors:", "\n", "\n")
        cat("Test Statistic =", round(statistic[1], digits=digits), "\n")
        cat("Degrees of Freedom =", df, "\n")
        cat("P-value of test =", round(p.value[1], digits=digits), "\n", "\n")
        cat("       Test for Equality of Dispersions:", "\n", "\n")
        cat("Test Statistic =", round(statistic[2], digits=digits), "\n")
        cat("Degrees of Freedom =", df, "\n")
        cat("P-value of test =", round(p.value[2], digits=digits), "\n", "\n")
    } else {
        cat("\n")
        cat("       Test for Equality of Polar Vectors:", "\n", "\n")
        cat("Test Statistic =", round(statistic[1], digits=digits), "\n")
        cat("Degrees of Freedom =", df, "\n")
        cat("Level", alpha, "critical value =", round(qchisq(1 - alpha, df), digits=digits), "\n")
        if (statistic[1] > qchisq(1 - alpha, df)) {
            cat("Reject null hypothesis of equal polar vectors", "\n", "\n")
        } else { 
                    cat("Do not reject null hypothesis of equal polar vectors", "\n", "\n")
                }
        cat("       Test for Equality of Dispersions:", "\n", "\n")
        cat("Test Statistic =", round(statistic[2], digits=digits), "\n")
        cat("Degrees of Freedom =", df, "\n")
        cat("Level", alpha, "critical value =", round(qchisq(1 - alpha, df), digits=digits), "\n")
        if (statistic[2] > qchisq(1 - alpha, df)) {
            cat("Reject null hypothesis of equal dispersions", "\n", "\n")
        } else {
                    cat("Do not reject null hypothesis of equal dispersions", "\n", "\n")
                }
    }
        invisible(x)
}
