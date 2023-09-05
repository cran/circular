.onAttach <- function(library, pkg)
{
  # Rv <- R.Version()
  # if(Rv$major < 2 |(Rv$major == 2 && Rv$minor < 2.0))
  #  stop("This package requires R 2.2.0 or later")
  if(interactive())
  {
    meta <- packageDescription("circular")
    # packageStartupMessage(
    #      "Package 'circular', ", meta$Version, " (", meta$Date, "). ",
    #      "Type 'help(Circular)' for summary information.")
  }
  invisible()
}
