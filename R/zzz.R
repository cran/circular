.First.lib <- function(lib, pkg) {
  library.dynam("circular", pkg, lib)
  cat("This is version 0.3 of circular package under tests \n")
  cat("Please report any bugs or comments to <Claudio Agostinelli> claudio@unive.it \n")
  cat("The package redefine how function 'density' and 'var' works \n")
  cat("In particular, for 'var' function (try 'methods(var)') \n notice that 'var.default' is an alias for the original 'var' function \n and that a method for data.frame is available.\n")
}
