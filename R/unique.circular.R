unique.circular <- function (x, incomparables = FALSE, ...) {
  if (!is.logical(incomparables) || incomparables) 
     .NotYetUsed("incomparables != FALSE")
  z <- .Internal(unique(x))
  circularp(z) <- circularp(x)
  class(z) <- class(x)
  return(z)
}
