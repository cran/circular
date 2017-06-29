regularize.values <- function (x, y, ties) {
    x <- xy.coords(x, y, setLab = FALSE)
    y <- x$y
    x <- x$x
    if (any(na <- is.na(x) | is.na(y))) {
        ok <- !na
        x <- x[ok]
        y <- y[ok]
    }
    nx <- length(x)
    if (!identical(ties, "ordered")) {
        o <- order(x)
        x <- x[o]
        y <- y[o]
        if (length(ux <- unique(x)) < nx) {
            if (missing(ties)) 
                warning("collapsing to unique 'x' values")
            y <- as.vector(tapply(y, match(x, x), ties))
            x <- ux
            stopifnot(length(y) == length(x))
        }
    }
    list(x = x, y = y)
}
