unique.circular <- function (x, incomparables = FALSE, fromLast = FALSE, ...) {
  if (!is.logical(incomparables) || incomparables) 
    .NotYetUsed("incomparables != FALSE")
    if (is.na(fromLast <- as.logical(fromLast[1]))) 
      stop("'fromLast' must be TRUE or FALSE")
  z <- .Internal(unique(x, fromLast))
  circularp(z) <- circularp(x)
  class(z) <- class(x)
  return(z)
}
