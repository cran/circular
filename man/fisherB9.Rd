\name{fisherB9}
\alias{fisherB9}
\alias{fisherB9direction}
\alias{fisherB9frequency}
\alias{fisherB9c}

\title{B.9 Dance directions of bees}

\usage{
data(fisherB9)
data(fisherB9c)
}

\description{
Dance directions of 279 honey bees viewing a zenith patch of artificially polarised light.
}
\format{
\code{fisherB9} a vector of 279 observations (in degrees).
\code{fisherB9c} contains the same observations in a circular objects.
}

\source{
Adapted by Prof. N.I. Fisher from 
R. Wehner & S. Strasser (1985) The POL area of the honey bee's eye: behavioural evidence. Physiol. Entomol. 10, 337-49. Pag. 346. 
}
\seealso{
N.I. Fisher (1993) Statistical analysis of circular data. Cambridge University Press. Pag. 244.
}

\examples{
data(fisherB9c)
plot(fisherB9c, stack=TRUE, shrink=1.5)
}

\keyword{datasets}
