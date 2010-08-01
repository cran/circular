\name{fisherB13}
\alias{fisherB13}
\alias{fisherB13c}

\title{B.13: Orientations of termite mounds}

\usage{
data(fisherB13)
data(fisherB13c)
}

\description{
Orientations of termite mounds of Amitermes laurensis at 14 sites in
Cape York Penisula, North Queensland. 
}
\format{
\code{fisherB13} a list of 14 datasets (axes in degrees) at several locations.
\code{fisherB13c} contains the same observations in a circular objects.}
\details{
  Set 1: n=100, Latitude -15'43'', Longitude 144'42''
  Set 2: n=50, Latitude -15'32'', Longitude 144'17''  
}
\source{
A.V. Spain, T. Okello-Oloya and R.D. John (1983) Orientation of the
termitaria of two species of Amitermes (Isoptera:Termitinae) from
Northern Queensland. Aust. J. Zoo. (31):167-177.
}
\references{
N.I. Fisher (1993) Statistical analysis of circular data. Cambridge University Press. Pag. 246.
}

\examples{
data(fisherB13c)
plot(fisherB13c$set1, stack=TRUE, shrink=1.5)
}

\keyword{datasets}
