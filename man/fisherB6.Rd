\name{fisherB6}
\alias{fisherB6}
\alias{fisherB6c}

\title{B.6 Cross-bed azimuths of palaeocurrents}

\usage{
data(fisherB6)
data(fisherB6c)
}

\description{
Set of cross-bed azimuths of palaeocurrents measured in the Belford Anticline (New South Wales).
}
\format{
\code{fisherB6} is a list (in degrees).
\code{fisherB6c} contains the same observations in a circular objects.
}

\source{
Fisher, N.I. & Powell C. McA. (1989) Statistical analysis of two-dimensional palaeocurrent data: Methods and examples. Aust. J. Earth Sci. 36, 91-107.
}
\seealso{
N.I. Fisher (1993) Statistical analysis of circular data. Cambridge University Press. Pag. 242.
}

\examples{
data(fisherB6c)
res <- plot(fisherB6c$set1)
points(fisherB6c$set2, col=2, plot.info=res)
points(fisherB6c$set3, col=3, plot.info=res)
}

\keyword{datasets}
