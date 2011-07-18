\name{wallraff.test}
\title{Wallraff Test of Angular Distances}
\alias{wallraff.test}
\alias{wallraff.test.default}
\alias{wallraff.test.list}
\alias{wallraff.test.formula}

\description{
Performs the Wallraff test of angular distances or angular dispersion around the mean.
}

\usage{
wallraff.test(x, ...)

\method{wallraff.test}{default}(x, group, ref=NULL, ...)

\method{wallraff.test}{list}(x, ref=NULL, ...)

\method{wallraff.test}{formula}(formula, data, ref=NULL, ...)
}

\arguments{
	\item{x}{a vector of angles (coerced to class \code{\link{circular}}) or a list of such angles. When \code{x} is a list, its elements are taken as the samples to be compared.}
	\item{group}{a vector or factor object giving the group for the corresponding elements of \code{x}. Ignored if \code{x} is a list}
	\item{formula}{a formula of the form \code{lhs ~ rhs} where \code{lhs} is a vector of angles and \code{rhs} a vector or factor giving the corresponding groups.}
	\item{data}{an optional data.frame containing the variables in the formula \code{\link{formula}}.}
	\item{ref}{a vector of angles used as reference to compute the angular distances from, in each group. It should contain as many elements as there are groups, \emph{in the same order}.
	
	If \code{x} is a list, the order is the order of the elements of the list.
	
	In the default or formula interfaces, if the grouping vector is a factor, the order is the order of its levels; if the grouping vector is not a factor, it is coerced as such but with levels in the order of their appearance in the original vector. In this case a warning is issued to make sure the order of \code{ref} is correct.
	
	If \code{ref} has less elements than the number of groups (typically one: a common reference for all groups), it is recycled to match the number of groups.
	
	If \code{ref} is \code{NULL} (the default), the mean angle of each group will be used as reference. In this situation, the Wallraff test becomes a comparison of angular dispersion around the mean.}
	\item{\dots}{further arguments passed to or from other methods.}
}

\details{
The Wallraff test of angular distances between two or more groups is performed and the results are printed. The null hypothesis is that distances are equal across groups.

The test proceeds by computing the angular distances from a reference angle, in each group. The angular distance between two angles is the circular range and is computed with \code{\link{range.circular}}. Then the distances are compared with a usual rank sum test (Kruskal-Wallis, \code{\link{kruskal.test}}). When there are only two groups, the Wilcoxon-Mann-Whitney test could be used but \code{\link{wilcox.test}} without continuity correction for the p-value is equivalent to \code{\link{kruskal.test}} so only \code{\link{kruskal.test}} is used here.

The Wallraff test is most frequently used to compare angular dispersion around the mean, between samples. In this case, the reference angle is the mean angle of each sample. This is the default here, when no reference angles are provided.

All angles should be of class \code{\link{circular}} and will be coerced as such with the default parameters if they are not. An exception are the reference angles in \code{ref}. For ease of use, those can be only numeric and are then considered to be in the same angular reference as \code{x}.
}

\value{
A list with class \code{"htest"} containing the following components:

	\item{statistic}{the chi-squared statistic from \code{\link{kruskal.test}}.}
	\item{parameter}{the degrees of freedom for the chi-squared statistic.}
	\item{p.value}{the p-value for the test.}
	\item{method}{a character string containing the name of the test.}
	\item{data.name}{a character string giving the name(s) of the data.}
}

\author{Jean-Olivier Irisson}

\references{
Batschelet, E (1981). Circular Statistics in Biology. chap. 6.10, p. 124

Zar, J H (2010). Biostatistical analysis. sec. 27.7-8, p. 643
}

\seealso{
\code{\link{kruskal.test}} for the Kruskal-Wallis rank sum test used on the angular distances.

\code{\link{wilcox.test}} for the two samples alternative to the Kruskal-Wallis test.
}

\examples{

# Homing of pigeons
# Example used in Batschelet (1981)
data <- list(
  control = circular(c(70, 80, 80, 85, 85, 90, 95, 95), 
        units="degrees", template="geographics"),
  experimental = circular(c(5, 5, 15, 55, 55, 65, 105, 120, 340), 
            units="degrees", template="geographics")
)

# compare the angular dispersion between the two groups
wallraff.test(data)

# compare the homing performance
# home azimuth is 40 degrees for both groups
wallraff.test(data, 
              ref = circular(c(40, 40), units="degrees",
                                 template="geographics")
             )
# we could have more simply used
wallraff.test(data, ref=40)
# because ref is automatically repeated and considered
# in the same circular reference as the data

}

\keyword{htest}
