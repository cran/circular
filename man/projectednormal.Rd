\name{projected bivariate normal on the circle}
\alias{dpnorm}
\alias{rpnorm}
\title{
Projected bivariate normal on the circle  
}
\description{
The projected normal distribution provides a flexible distribution for circular data, e.g., asymmetry and possible bimodality.
}
\usage{
dpnorm(x, mu, sigma, log = FALSE)
rpnorm(n, mu, sigma, control.circular=list())
}
\arguments{
  \item{x}{a vector. The \code{x} and \code{q} objects are coerced to class \code{\link{circular}}.
}
  \item{n}{number of observations.} 
  \item{mu}{the mean vector of the bivariate normal.}
  \item{sigma}{the 2x2 variance and covariance matrix of the bivariate normal.}
  \item{log}{logical. If \code{TRUE} the log of the density is reported.}
  \item{control.circular}{the attribute of the resulting object.}  
}

\value{\code{dpnorm} gives the density, \code{rpnorm} generates random deviates.}
\references{
S.R. Jammalamadaka and A. SenGupta (2001). Topics in Circular Statistics, Section 2.2.4, World Scientific Press, Singapore.
K.V. Mardia (1972). Statistics of Directional Data. Academic Press. London and New York.    
F. Wang and A.E. Gelfand (2013). Directional data analysis under the general projected normal distribution. Stat Methodol. 10(1):113-127. doi:10.1016/j.stamet.2012.07.005.
}
\author{
Claudio Agostinelli
}

\examples{
data1 <- rpnorm(100, mu=c(0,0), sigma=diag(2),
  control.circular=list(units="degrees")) # Uniform on the circle
plot(data1)

ff <- function(x) dpnorm(x, mu=c(0,0), sigma=diag(2)) # Uniform on the circle
curve.circular(ff, join=TRUE,
  main="Density of a Projected Normal Distribution \n mu=(0,0), sigma=diag(2)")

ff <- function(x) dpnorm(x, mu=c(1,1), sigma=diag(2)) # Unimodal
curve.circular(ff, join=TRUE, xlim=c(-1, 2.3),
  main="Density of a Projected Normal Distribution \n mu=(1,1), sigma=diag(2)")

sigma <- matrix(c(1,0.9,0.9,1), nrow=2)  
ff <- function(x) dpnorm(x, mu=c(0.5,0.5), sigma=sigma) # Bimodal
curve.circular(ff, join=TRUE, xlim=c(-1, 2.3),
  main="Density of a Projected Normal Distribution \n mu=(0.5,0.5), rho=0.9")
}

\keyword{distribution}
