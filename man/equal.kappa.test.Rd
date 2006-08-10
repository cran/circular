\name{equal.kappa.test}
\alias{equal.kappa.test}
\alias{print.equal.kappa.test}
\title{Equal Kappa Test}
\description{
This function tests for the homogeneity of concentration parameters for multiple 
samples of directional data.
}
\usage{
equal.kappa.test(x, group)
\method{print}{equal.kappa.test}(x, digits = max(3, getOption("digits") - 3), ...)

}

\arguments{
  \item{x}{a vector of class \code{circular}.}
  \item{group}{a vector identifying the groups or samples.}
  \item{digits}{the number of digits to be printed.}
  \item{\dots}{additional arguments.}
}

\value{
  An object of class \code{equal.kappa.test} with the following
  components:
  
  \item{kappa}{concentration parameter for each sample.}
  \item{kappa.all}{concentration parameter of all samples combined.}
  \item{rho}{mean resultant length for each sample.}
  \item{rho.all}{mean resultant length of all samples combined.}
  \item{df}{degrees of freedom for chi-squared distribution.}
  \item{statistic}{the value of the chi-squared test statistic.}
  \item{p.value}{the p.value of the test statistic.}
  \item{call}{match.call()}
}

\details{
The samples are assumed to have been drawn from von Mises populations. 
The null hypothesis tested is that all populations sampled have the same 
concentration parameter, kappa. 

When the pooled data has high concentration, sample mean resultant length 
above 0.70, Bartlett's test is used. For less concentrated pooled data, 
variance-stabilizing transformations are used to improve normal 
approximations needed to arrive at an approximate chi-squared 
test statistic (see references below). For pooled sample mean resultant 
length below 0.45, it is possible that individually a sample may in fact have 
quite a large sample mean resultant length. In this case, it is possible that 
the variance-stabilizing transformation involving the inverse sine function 
is passed a value outside of -1,1. If this occurs, the function will 
automatically use Bartlett's test and issue a warning to that effect. 
}

\references{
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, 
Section 5.3, World Scientific Press, Singapore. 

Mardia, K. and Jupp, P. (1999). Directional Statistics, Section 7.4, John 
Wiley and Sons, England. 
}

\author{Claudio Agostinelli and Ulric Lund}

\examples{

x <- c(rvonmises(50, circular(0), 1), rvonmises(100, circular(pi/3), 10))
group <- c(rep(0, 50), rep(1, 100)) 

equal.kappa.test(x, group)

}

\keyword{htest}
