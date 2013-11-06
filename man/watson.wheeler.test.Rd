\name{watson.wheeler.test}
\title{Watson-Williams Test of Homogeneity of Means}
\alias{watson.wheeler.test}
\alias{watson.wheeler.test.default}
\alias{watson.wheeler.test.list}
\alias{watson.wheeler.test.formula}

\description{
Performs the Watson-Wheeler test for homogeneity on two or more samples of circular data. 
}

\usage{
watson.wheeler.test(x, ...)

\method{watson.wheeler.test}{default}(x, group, ...)

\method{watson.wheeler.test}{list}(x, ...)

\method{watson.wheeler.test}{formula}(formula, data, ...)
}

\arguments{
	\item{x}{a vector of angles (coerced to class \code{\link{circular}}) or a list of such angles.}
	\item{group}{a vector or factor object giving the groups for the corresponding elements of \code{x}. Ignored if \code{x} is a list}
	\item{formula}{a formula of the form \code{lhs ~ rhs} where \code{lhs} is a vector of angles and \code{rhs} a vector or factor giving the corresponding groups.}
	\item{data}{an optional data.frame containing the variables in the formula \code{\link{formula}}.}
	\item{\dots}{further arguments passed to or from other methods.}
}

\value{
A list with class \code{"htest"} containing the following components:

	\item{statistic}{W, the statistic of the test, which is approximately distributed as a chi-squared.}
	\item{parameter}{the degrees of freedom for the chi-squared approximation of the statistic.}
	\item{p.value}{the p-value for the test.}
	\item{method}{a character string containing the name of the test.}
	\item{data.name}{a character string giving the name(s) of the data.}
}

\details{
The Watson-Wheeler (or Mardia-Watson-Wheeler, or uniform score) test is a non-parametric test to compare two or several samples. The difference between the samples can be in either the mean or the variance.

The \emph{p}-value is estimated by assuming that the test statistic follows a chi-squared distribution. For this approximation to be valid, all groups must have at least 10 elements.

In the default method, \code{x} is a vector of angles and \code{group} must be a vector or factor object of the same length as \code{x} giving the group for the corresponding elements of \code{x}.

If \code{x} is a list, its elements are taken as the samples to be compared.

In the \code{\link{formula}} method, the angles and grouping elements are identified as the left and right hand side of the formula respectively.

All angles should be of class \code{\link{circular}} and will be coerced as such if they are not.
}

\author{Jean-Olivier Irisson}

\references{
Batschelet, E (1981). Circular Statistics in Biology. chap 6.3, p. 104

Zar, J H (1999). Biostatistical analysis. section 27.5, p. 640
}

\examples{

# Example used in Zar (1999)
x1 <- circular(c(35, 45, 50, 55, 60, 70, 85, 95, 105, 120),
  units="degrees", template="geographics")
x2 <- circular(c(75, 80, 90, 100, 110, 130, 135, 140, 150, 160, 165),
  units="degrees", template="geographics")

watson.wheeler.test(list(x1,x2))
}

\keyword{htest}
