\name{fisherB10}
\alias{fisherB10}
\alias{fisherB10c}

\title{B.10 Directions of desert ants}

\usage{
data(fisherB10)
data(fisherB10c)
}

\description{
Directions of 11 long-legged desert ants (Cataglyphis fortis) after one eye on each ant was 'trained' to learn the ant's home direction, then covered and the other eye uncovered.
}
\format{
\code{fisherB10} is a list (in degrees).
\code{fisherB10c} contains the same observations in a circular objects.
}

\source{
Personal communication of Prof. Dr. R. Wehner to Prof. N.I. Fisher, experiment described in 

R. Wehner and M. Muller (1985) Does interocular transfer occur in visual navigation by ants? Nature, 315, 228-9.
}
\seealso{
N.I. Fisher (1993) Statistical analysis of circular data. Cambridge University Press. Pag. 244-245.
}

\examples{
data(fisherB10c)
res <- plot(fisherB10c$set1)
points(fisherB10c$set2, col=2, plot.info=res)
points(fisherB10c$set3, col=3, plot.info=res)
}

\keyword{datasets}
