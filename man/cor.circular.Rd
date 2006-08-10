\name{cor.circular}
\title{Correlation Coefficient for Angular Variables}
\alias{cor.circular}

\description{
Computes a circular version of the Pearson's product moment correlation, and performs a significance test if requested. 
}

\usage{
cor.circular(x, y=NULL, test=FALSE)
}

\arguments{
\item{x}{vector or matrix of circular data.}
\item{y}{vector or matrix of circular data.}
\item{test}{if \code{test == TRUE}, then a significance test for the correlation coefficient is computed.}
}

\value{
Returns a vector or a matrix of a circular version of the Pearson's
product moment correlation, if \code{test == TRUE} then a list is
reported with statistic and p.value, the test statistic and p-value respectively, for testing significance of the correlation coefficient.
}

\details{
The correlation coefficient is computed like Pearson's product moment correlation for two linear variables X and Y.  In the computational formula,  however, (xi - xbar) and (yi - ybar) are replaced by sin(xi - xbar) and sin(yi - ybar), where xbar and ybar in the second two expressions are the mean directions of the samples.
}

\author{Claudio Agostinelli and Ulric Lund}

\references{
Jammalamadaka, S. Rao and SenGupta, A. (2001). Topics in Circular Statistics, Section 8.2, World Scientific Press, Singapore.

Jammalamadaka, S. and Sarma, Y. (1988).  A correlation coefficient for angular variables.  Statistical Theory and Data Analysis 2. North Holland:  New York.
}

\examples{
# Generate two circular data sets, and compute their correlation.
x <- rvonmises(n=50, mu=circular(0), kappa=3)
y <- x + rvonmises(n=50, mu=circular(pi), kappa=10)
cor.circular(x, y, test=TRUE)
}

\keyword{univar}
