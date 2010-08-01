\name{fisherB7}
\alias{fisherB7}
\alias{fisherB7c}

\title{B.7 Movements of ants}

\usage{
data(fisherB7)
data(fisherB7c)
}

\description{
Directions chisen by 100 ants in response to an evenly illuminated black targed placed as shown.
}
\format{
\code{fisherB7} a vector of 100 observations (in degrees).
\code{fisherB7c} contains the same observations in a circular objects.
}

\source{
Randomly selected values from 
Jander, R. (1957) Die optische Richtangsorientierung der roten Waldameise (Formica rufa. L.) Z. vergl. Physiologie 40, 162-238. Figure 18A. 
}
\seealso{
N.I. Fisher (1993) Statistical analysis of circular data. Cambridge University Press. Pag. 243.
}

\examples{
data(fisherB7c)
plot(fisherB7c, zero=pi/2, rotation='clock', stack=TRUE)
}

\keyword{datasets}
