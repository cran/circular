\name{fisherB11}
\alias{fisherB11}
\alias{fisherB11c}

\title{B.11 Movements of sea stars}

\usage{
data(fisherB11)
data(fisherB11c)
}

\description{
Resultant directions of 22 sea stars 11 days after being displaced from their natural habitat.
}
\format{
\code{fisherB11} a vector of 22 observations (in degrees).
\code{fisherB11c} contains the same observations in a circular objects.
}

\source{
G.J.G. Upton and B. Fingleton (1989) Spatial Data Analysis by Example. Volume 2. Categorical and Directional Data. New York: John Wiley
as adapted from
B. Pabst and H. Vicentini (1978) Dislocation experiments in the migrating seastar. Astropecten jonstoni. Marine Biology 48, 271-8. 
}
\seealso{
N.I. Fisher (1993) Statistical analysis of circular data. Cambridge University Press. Pag. 245.
}

\examples{
data(fisherB11c)
plot(fisherB11c, stack=TRUE, shrink=1.5)
}

\keyword{datasets}
