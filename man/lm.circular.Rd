\name{lm.circular}    
\title{Circular-Circular Regression}
\alias{lm.circular}

\description{
Fits a regression model for a circular dependent and circular independent variable.  
}

\usage{
lm.circular(y, x, order=1, level=0.05)
}

\arguments{
\item{y}{vector of data for the dependent circular variable.}
\item{x}{vector of data for the independent circular variable.}
\item{order}{order of trigonometric polynomial to be fit. Order must be an integer value.  By default, order=1.}
\item{level}{level of the test for the significance of higher order trigonometric terms.}
}

\value{
A list is returned with the following components: 

\item{call}{match.call().}
\item{rho}{square root of the average of the squares of the estimated conditional concentration parameters of y given x.}
\item{fitted}{fitted values of the model.}
\item{data}{matrix whose columns correspond to x and y.}
\item{residuals}{circular residuals of the model.}
\item{coefficients}{matrix whose entries are the estimated coefficients of the model. The first column corresponds to the coefficients of the model predicting the cosine of y, while the second column contains the estimates for the model predicting the sine of y.  The rows of the matrix correspond to the coefficients according to increasing trigonometric order.}
\item{p.values}{p-values testing whether the (order + 1) trigonometric terms are significantly different from zero.}  
\item{A.k}{is mean of the cosines of the circular residuals.}
\item{kappa}{assuming the circular residuals come from a von Mises distribution, kappa is the MLE of the concentration parameter.}
}

\details{
A trigonometric polynomial of x is fit against the cosine and sine of y.  The order of trigonometric polynomial is specified by order.  Fitted values of y are obtained by taking the inverse tangent of the predicted values of sin(y) devided by the predicted values of cos(y).  Details of the regression model can be found in Sarma and Jammalamadaka (1993).
}

\author{Claudio Agostinelli and Ulric Lund}

\references{
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 8.3, World Scientific Press, Singapore.

Sarma, Y. and Jammalamadaka, S. (1993).  Circular Regression.  Statistical Science and Data Analysis, 109–128.  Proceeding of the Thrid Pacific Area Statistical Conference.  VSP: Utrecht, Netherlands.
}

\examples{
# Generate a data set of dependent circular variables.
x <- circular(runif(50, 0, 2*pi))
y <- atan(0.15*cos(x) + 0.25*sin(x), 0.35*sin(x)) + rvonmises(n=50, mu=0, kappa=5)

# Fit a circular regression model.
circ.lm <- lm.circular(y, x, order=1)
# Obtain a crude plot a data and fitted regression line.
plot.default(x, y)
circ.lm$fitted[circ.lm$fitted>pi] <- circ.lm$fitted[circ.lm$fitted>pi] - 2*pi 

points.default(x[order(x)], circ.lm$fitted[order(x)], type='l')
}

\keyword{models}
