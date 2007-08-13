\name{aov.circular}
\alias{aov.circular}
\alias{print.aov.circular}
\title{Analysis of Variance for circular data}
\description{

}

\usage{
aov.circular(x, group, kappa = NULL,
    method = c("F.test", "LRT"), F.mod = TRUE, control.circular=list())
\method{print}{aov.circular}(x, digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{x}{a vector of class \code{circular}.}
  \item{group}{a vector identifying the groups or samples.}
  \item{kappa}{the common value of the concentration parameter. Used
    when \code{method} is \code{"LRT"}.  If left unspecified (by default) 
    the maximum likelihood estimate of kappa is computed and used in the test statistic.}
  \item{method}{the test statistic to use; either a high-concentration F-test or a 
    likelihood ratio test.}
  \item{F.mod}{logical; if \code{TRUE}, the AOV F-statistic is modified by a factor of 1+3/8k 
    to improve its approximate F distribution. Default is \code{TRUE}.}
  \item{control.circular}{the coordinate system used in the output for the objects \code{mu} and \code{mu.all}. See \code{\link{circular}} for details.}
  \item{digits}{the number of digits to be printed.}
  \item{\dots}{additional arguments.}
}

\value{
  An object of class \code{aov.circular} with the folowing components:

  \item{mu}{mean direction for each sample with class \code{circular}.}
  \item{mu.all}{mean direction of all samples combined with class \code{circular}.}
  \item{kappa}{concentration parameter for each sample.}
  \item{kappa.all}{concentration parameter for all samples combined.}
  \item{rho}{mean resultant length for each sample.}
  \item{rho.all}{mean resultant length for all samples combined.}
  \item{method}{the test statistic used.}
  \item{df}{degrees of freedom.}
  \item{statistic}{the value of the test statistic.}
  \item{p.value}{the p.value of the test statistic.}
  \item{call}{match.call().}

  If the \code{method} is \code{"F.test"} then the object contains also:

  \item{SSE}{Sum of squares used in F-test.}
  \item{MSE}{Mean squares used in F-test.}
}

\details{
The samples are assumed to have been drawn from von Mises populations with 
equal concentration parameter, kappa. The null hypothesis being tested is 
that all populations also have the same mean direction. 

If \code{method} is \code{"F.test"} a high concentration F-test makes use 
of a decomposition of total sample variation into between groups and within 
groups variation, analogous to the one-way classification analysis of variance 
for linear data.  Stephens (1972) presented an improved modification to the 
F-test derived from this decomposition.  This is implemented when \code{F.mod} 
is \code{TRUE}.

A likelihood ratio test is performed when \code{method} is \code{"LRT"}. This 
function uses the test statistic presented by Cordeiro, Paula, and Botter (1994) which has 
an approximate chi-squared distribution. If the common concentration parameter 
is known, it can be specified and used in the computation of the test statistic. 
Otherwise, the maximum likelihood estimate of the common concentration parameter is used.
}

\references{
Cordeiro, G., Paula, G. and Botter, D. (1994). Improved likelihood ratio tests for 
dispersion models. International Statistical Review, 62, 257-274. 

Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 5.3, 
World Scientific Press, Singapore. 

Mardia, K. and Jupp, P. (1999). Directional Statistics, Section 7.4, John 
Wiley and Sons, England. 

Stephens, M. (1972). Multi-sample tests for the von Mises distribution. Technical 
Report 190, Department of Statistics, Stanford University. 
}

\author{Claudio Agostinelli and Ulric Lund}

\examples{
x <- c(rvonmises(50, circular(0), 1), rvonmises(100, circular(pi/3), 10))
group <- c(rep(0, 50), rep(1, 100))

aov.circular(x, group)
aov.circular(x, group, method="LRT")
}
\keyword{models}
