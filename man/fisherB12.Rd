\name{fisherB12}
\alias{fisherB12}
\alias{fisherB12c}

\title{B.12: Vanishing directions of homing pigeons}

\usage{
data(fisherB12)
data(fisherB12c)
}

\description{
Vanishing directions of 15 homing pigeons, released just over
16 kilometres Northwest of their loft.
}
\format{
\code{fisherB12} a vector of 15 observations (in degrees).
\code{fisherB12c} contains the same observations in a circular objects.}

\source{
Schmidt-Koenig (1963) On the role of the loft, the distance and site of
release in pigeon homing (the "cross-loft experiment"). Biol. Bull. (125)154-164.
}
\references{
N.I. Fisher (1993) Statistical analysis of circular data. Cambridge University Press. Pag. 245.
}

\examples{
data(fisherB12c)
plot(fisherB12c, stack=TRUE, shrink=1.5)
}

\keyword{datasets}
