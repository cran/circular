\name{watson.williams.test}
\title{Watson-Williams Test of Homogeneity of Means}
\alias{watson.williams.test}
\alias{watson.williams.test.default}
\alias{watson.williams.test.list}
\alias{watson.williams.test.formula}

\description{
Performs the Watson-Williams test for homogeneity of means between several samples of circular data. 
}

\usage{
watson.williams.test(x, ...)

\method{watson.williams.test}{default}(x, group, ...)

\method{watson.williams.test}{list}(x, ...)

\method{watson.williams.test}{formula}(formula, data, ...)
}

\arguments{
	\item{x}{a vector of angles (coerced to class \code{\link{circular}}) or a list of such angles.}
	\item{group}{a vector or factor object giving the group for the corresponding elements of \code{x}. Ignored if \code{x} is a list}
	\item{formula}{a formula of the form \code{lhs ~ rhs} where \code{lhs} is a vector of angles and \code{rhs} a vector or factor giving the corresponding groups.}
	\item{data}{an optional data.frame containing the variables in the formula \code{\link{formula}}.}
	\item{\dots}{further arguments passed to or from other methods.}
}

\value{
A list with class \code{"htest"} containing the following components:

	\item{statistic}{the F statistic of the test.}
	\item{parameter}{the degrees of freedom for the F statistic.}
	\item{p.value}{the p-value for the test.}
	\item{estimate}{a vector of the means of each group.}
	\item{method}{a character string containing the name of the test.}
	\item{data.name}{a character string giving the name(s) of the data.}
}

\details{
The Watson-Williams test for the homogeneity of means between two or more groups is performed and the results are printed. The null hypothesis is that means are equal across groups.

The assumptions are that: (1) the samples are drawn from populations with a von Mises distribution; (2) the parameter of concentration has the same value in all populations; (3) this parameter is sufficiently large (i.e. > 1). Assumptions 2 and 3 are checked and a warning is issued if they are not met.

In the default method, \code{x} is a vector of angles and \code{group} must be a vector or factor object of the same length as \code{x} giving the group for the corresponding elements of \code{x}.

If \code{x} is a list, its elements are taken as the samples to be compared.

In the \code{\link{formula}} method, the angles and grouping elements are identified as the left and right hand side of the formula respectively.

All angles should be of class \code{\link{circular}} and will be coerced as such if they are not.
}

\author{Jean-Olivier Irisson}

\references{
Batschelet, E (1981). Circular Statistics in Biology. chap. 6.2, p. 99

Mardia, KV and Jupp, PE (2000). Directional statistics. p. 135
}

\examples{

# Ant orientation from Duelli and Wehner (1973)
# Example used in Batschelet (1981)
data <- list(
  exp = circular(rep(c(-20, -10, 0), c(1,7,2)), 
        units="degrees", template="geographics"),
  control = circular(rep(c(-10, 0, 10, 20), c(3,3,3,1)), 
            units="degrees", template="geographics")
)

watson.williams.test(data)

}

\keyword{htest}
