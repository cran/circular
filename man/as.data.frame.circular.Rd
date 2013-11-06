\name{as.data.frame.circular}
\title{as.data.frame.circular}
\alias{as.data.frame.circular}

\description{
  This function is a method of \code{\link{as.data.frame}} for a circular object.
}
\usage{
\method{as.data.frame}{circular}(x, row.names = NULL, optional = FALSE, ...)
}
\arguments{
  \item{x}{object of class \code{\link{circular}}.}
  \item{row.names}{\code{NULL} or a character vector giving the row
    names for the data frame. Missing values are not allowed.}
  \item{optional}{logical; if \code{TRUE} setting row names is
      optional.}
  \item{\dots}{additional arguments to be passed to or from methods.}
}

\author{Claudio Agostinelli}

\keyword{classes}
\keyword{methods}
