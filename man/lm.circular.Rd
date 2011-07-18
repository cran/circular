\name{lm.circular}    
\title{Circular-Circular and Circular-Linear Regression}
\alias{lm.circular}
\alias{lm.circular.cc}
\alias{lm.circular.cl}
\alias{print.lm.circular.cl}

\description{
Fits a regression model for a circular dependent and circular
independent variable or for a circular dependent and linear independent
variables.  
}

\usage{
lm.circular(..., type=c("c-c", "c-l"))
lm.circular.cc(y, x, order = 1, level = 0.05, control.circular = list())
lm.circular.cl(y, x, init = NULL, verbose = FALSE, tol = 1e-10, 
  control.circular = list())
\method{print}{lm.circular.cl}(x, digits = max(3, getOption("digits") - 3), 
  signif.stars= getOption("show.signif.stars"), ...)
}

\arguments{
\item{\dots}{arguments passed to \code{lm.circular.cc} or to
  \code{lm.circular.cl} depending on the value of \code{type}.}
\item{type}{if \code{type=="c-c"} then \code{lm.circular.cc} is called
  otherwise \code{lm.circular.cl} is called.}
\item{y}{vector of data for the dependent circular variable.}
\item{x}{vector of data for the independent circular variable if
  \code{type="c-c"} or \code{lm.circular.cc} is used otherwise a matrix
  or a vector containing the independent linear variables.}
\item{order}{order of trigonometric polynomial to be fit. Order must be
  an integer value.  By default, order=1. Used if \code{type="c-c"}.}
\item{level}{level of the test for the significance of higher order
  trigonometric terms. Used if \code{type="c-c"}.}
\item{control.circular}{the attribute of the resulting objects (\code{fitted}, \code{residuals} components in the case of \code{type=="c-c"} and \code{mu} and \code{se.mu}) otherwise.}
\item{init}{a vector with initial values of length equal to the columns
  of \code{x}.}
\item{verbose}{logical: if \code{TRUE} messages are printed while the
  function is running.}
\item{tol}{the absolute accuracy to be used to achieve convergence of the
  algorithm.}
\item{digits}{the number of digits to be printed.}
\item{signif.stars}{logical; if \code{TRUE}, P-values are additionally encoded
   visually as ``significance stars'' in order to help scanning of long
   coefficient tables. It defaults to the \code{show.signif.stars} slot of \code{\link{options}}.}
}

\value{

If \code{type=="c-c"} or \code{lm.circular.cc} is called directly an
object of class \code{lm.circular.cc} is returned with the following components:

\item{call}{the \code{\link[base]{match.call}} result.}
\item{rho}{square root of the average of the squares of the estimated
conditional concentration parameters of y given x.}
\item{fitted}{fitted values of the model of class \code{circular}.}
\item{data}{matrix whose columns correspond to x and y.}
\item{residuals}{circular residuals of the model of class \code{circular}.}
\item{coefficients}{matrix whose entries are the estimated coefficients of the
model. The first column corresponds to the coefficients of the model predicting
the cosine of y, while the second column contains the estimates for the model
predicting the sine of y.  The rows of the matrix correspond to the coefficients
according to increasing trigonometric order.}
\item{p.values}{p-values testing whether the (order + 1) trigonometric terms are
significantly different from zero.}  
\item{A.k}{is mean of the cosines of the circular residuals.}
\item{kappa}{assuming the circular residuals come from a von Mises
  distribution, kappa is the MLE of the concentration parameter.}

If \code{type=="c-l"} or \code{lm.circular.cl} is called directly an
object of class \code{lm.circular.cc} is returned with the following components:

\item{call}{the \code{\link[base]{match.call}} result.}
\item{x}{the independent variables.}
\item{y}{the dependent variable.}
\item{mu}{the circular mean of the dependent variable of class \code{circular}.}
\item{se.mu}{an estimated standard error of the circular mean with the same units of measure used for \code{mu}.} 
\item{kappa}{the concentration parameter for the dependent variable.}
\item{se.kappa}{an estimated standard error of the concentration parameter.}
\item{coefficients}{the estimated coefficients.}
\item{cov.coef}{covariance matrix of the estimated coefficients.}
\item{se.coef}{standard errord of the estimated coefficients.}
\item{log.lik}{log-likehood.}
\item{t.values}{values of the t statistics for the coefficients.}
\item{p.values}{p-values of the t statistics. Approximated values using
  Normal distribution.}
}

\details{
  If \code{type=="c-c"} or \code{lm.circular.cc} is called directly a
  trigonometric polynomial of x is fit against the cosine and sine of y.
  The order of trigonometric polynomial is specified by order.  Fitted
  values of y are obtained by taking the inverse tangent of the predicted
  values of sin(y) devided by the predicted values of cos(y).  Details of
  the regression model can be found in Sarma and Jammalamadaka (1993).
  
  If \code{type=="c-l"} or \code{lm.circular.cl} is called directly, 
  this function implements the homoscedastic version of the maximum 
  likelihood regression model proposed by Fisher and Lee (1992). The 
  model assumes that a circular response variable theta has a von Mises 
  distribution with concentration parameter kappa, and mean direction 
  related to a vector of linear predictor variables according to the 
  relationship: mu + 2*atan(beta'*x), where mu and beta are unknown 
  parameters, beta being a vector of regression coefficients. The 
  function uses Green's (1984) iteratively reweighted least squares 
  algorithm to perform the maximum likelihood estimation of kappa, mu, 
  and beta. Standard errors of the estimates of kappa, mu, and beta are 
  estimated via large-sample asymptotic variances using the information 
  matrix. An estimated circular standard error of the estimate of mu is 
  then obtained according to Fisher and Lewis (1983, Example 1). 
}

\author{Claudio Agostinelli and Ulric Lund}

\references{
Fisher, N. and Lee, A. (1992). Regression models for an angular response. 
Biometrics, 48, 665-677. 

Fisher, N. and Lewis, T. (1983). Estimating the common mean direction of 
several circular or spherical distributions with different dispersions. 
Biometrika, 70, 333-341. 

Green, P. (1984). Iteratively reweighted least squares for maximum 
likelihood estimation, and some robust and resistant alternatives. 
Journal of the Royal Statistical Society, B, 46, 149-192. 

Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics,
Section 8.3, World Scientific Press, Singapore.

Sarma, Y. and Jammalamadaka, S. (1993).  Circular Regression.  Statistical
Science and Data Analysis, 109-128.  Proceeding of the Thrid Pacific Area
Statistical Conference.  VSP: Utrecht, Netherlands.
}

\examples{
# Generate a data set of dependent circular variables.
x <- circular(runif(50, 0, 2*pi))
y <- atan2(0.15*cos(x) + 0.25*sin(x), 0.35*sin(x)) + 
  rvonmises(n=50, mu=circular(0), kappa=5)

# Fit a circular-circular regression model.
circ.lm <- lm.circular(y, x, order=1)
# Obtain a crude plot of the data and fitted regression line.
plot.default(x, y)
circ.lm$fitted[circ.lm$fitted>pi] <- circ.lm$fitted[circ.lm$fitted>pi] - 2*pi 

points.default(x[order(x)], circ.lm$fitted[order(x)], type='l')


# Fit a circular-linear regression model.
set.seed(1234)
x <- cbind(rnorm(10), rep(1, 10))
y <- circular(2*atan(c(x\%*\%c(5,1))))+rvonmises(10, mu=circular(0), kappa=100)
lm.circular(y=y, x=x, init=c(5,1), type='c-l', verbose=TRUE)
}

\keyword{models}
\keyword{regression}
