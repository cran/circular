
###############################################################
#                                                             #
#       Original Splus: Ulric Lund                            #
#       E-mail: ulund@calpoly.edu                             #
#                                                             #
###############################################################

#############################################################
#                                                           #
#   lm.circular function                                    #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: August, 01, 2003                                  #
#   Version: 0.1-1                                          #
#                                                           #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################

lm.circular <- function(y, x, order = 1, level = 0.05) {
    x <- as.circular(x)
    x <- conversion.circular(x, units="radians")
    attr(x, "circularp") <- attr(x, "class") <- NULL
    y <- as.circular(y)
    y <- conversion.circular(y, units="radians")
    attr(y, "circularp") <- attr(y, "class") <- NULL

    n <- length(x)
    cy <- cos(y)
    sy <- sin(y)
    order.matrix <- t(matrix(rep(c(1:order), n), ncol = n))
    cos.x <- cos(x * order.matrix)
    sin.x <- sin(x * order.matrix)
    cos.lm <- lm(cy ~ cos.x + sin.x)
    sin.lm <- lm(sy ~ cos.x + sin.x)
    cos.fit <- cos.lm$fitted
    sin.fit <- sin.lm$fitted
    g1.sq <- t(cos.fit) %*% cos.fit
    g2.sq <- t(sin.fit) %*% sin.fit
    rho <- sqrt((g1.sq + g2.sq)/n)
    y.fitted <- atan(sin.fit, cos.fit)
    Y1 <- cy
    Y2 <- sy
    ones <- matrix(1, n, 1)
    X <- cbind(ones, cos.x, sin.x)
    W <- cbind(cos((order + 1) * x), sin((order + 1) * x))
    M <- X %*% solve(t(X) %*% X) %*% t(X)
    I <- diag(n)
    H <- t(W) %*% (I - M) %*% W
    N <- W %*% solve(H) %*% t(W)
    cc <- n - (2 * order + 1)
    N1 <- t(Y1) %*% (I - M) %*% N %*% (I - M) %*% Y1
    D1 <- t(Y1) %*% (I - M) %*% Y1
    T1 <- cc * (N1/D1)
    N2 <- t(Y2) %*% (I - M) %*% N %*% (I - M) %*% Y2
    D2 <- t(Y2) %*% (I - M) %*% Y2
    T2 <- cc * (N2/D2)
    p1 <- 1 - pchisq(T1, 2)
    p2 <- 1 - pchisq(T2, 2)
    pvalues <- cbind(p1, p2)
    circ.lm <- list()
    circ.lm$call <- match.call()
    circ.lm$rho <- rho
    circ.lm$fitted <- y.fitted %% (2 * pi)
    circ.lm$x <- cbind(x, y)
    circ.lm$residuals <- (y - y.fitted) %% (2 * pi)
    circ.lm$coefficients <- cbind(cos.lm$coefficients, sin.lm$coefficients)
    circ.lm$p.values <- pvalues
    circ.lm$A.k <- mean(cos(circ.lm$residuals))
    circ.lm$kappa <- A1inv(circ.lm$A.k)
    if (pvalues[1] > level & pvalues[2] > level)
        circ.lm$message <- paste("Higher order terms are not significant at the ", level, " level", sep = "")
    else circ.lm$message <- paste("Higher order terms are significant at the ", level, " level", sep = "")
    return(circ.lm)
}
