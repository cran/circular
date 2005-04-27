\name{plot.lsfit.circle}
\alias{plot.lsfit.circle}

\title{Plot method for lsfit.circle function}
\description{
This is a plot method for objects of class \code{lsfit.circle}.
}

\usage{
\method{plot}{lsfit.circle}(x, add = FALSE, main = NULL, xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL, uin, tol = 0.04, plus.cex = 1, ...)
}

\arguments{
  \item{x}{an object of class \code{lsfit.circle}.}
  \item{add}{logical: if \code{TRUE} the plot is superimposed on the
    active device.}
  \item{main}{a main title for the plot.}
  \item{xlim}{the x limits (min,max) of the plot.}
  \item{ylim}{the y limits of the plot.}
  \item{xlab}{a label for the x axis.}
  \item{ylab}{a label for the x axis.}
  \item{uin}{desired values for the units per inch parameter. If of length 1, the desired units per inch on the x axis.}
  \item{tol}{proportion of white space at the margins of plot.}
  \item{plus.cex}{dimension of the cross in the center of the circle.}
  \item{\dots}{further arguments passed to the next method.}
}

\author{Claudio Agostinelli and Ulric Lund}

\seealso{\code{\link{lsfit.circle}}}
\examples{

data(coope)
res <- lsfit.circle(x=x.coope, y=y.coope)
plot(res)
}

\keyword{hplot}
