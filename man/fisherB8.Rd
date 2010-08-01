\name{fisherB8}
\alias{fisherB8}
\alias{fisherB8c}

\title{B.8 Orientations of pebbles}

\usage{
data(fisherB8)
data(fisherB8c)
}

\description{
Horizontal axes of 100 outwash pebbles fromo a late Wisconsin outwash terrace alogn Fox river, near Cary, Illinois
}
\format{
\code{fisherB8} a vector of 100 observations (in degrees).
\code{fisherB8c} contains the same observations in a circular objects.
}

\source{
Mardia, K.V. (1972) Statistics of Directional Data. London: Academic Press. Table 1.6
adapted from
Krumbein W.C. (1939) Preferred orientations of pebbles in sedimentary deposits. J. Geol. 47, 673-706.
}
\seealso{
N.I. Fisher (1993) Statistical analysis of circular data. Cambridge University Press. Pag. 243.
}

\examples{
data(fisherB8c)
plot(fisherB8c, stack=TRUE, shrink=1.5)
}

\keyword{datasets}
