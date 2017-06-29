## .First.lib <- function(lib, pkg) {
##   library.dynam("circular", pkg, lib)
## }

.onLoad <- function(libname, pkgname) {
  data("rao.table", package=pkgname, envir=parent.env(environment()))
}  