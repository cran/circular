\name{fisherB1}
\alias{fisherB1}
\alias{fisherB1c}

\title{B.1 Arrival times at an intensive care unit}

\usage{
data(fisherB1)
data(fisherB1c)
}

\description{
Arrival time on a 24-hour clock of 254 patients at an intensive care unit, over a period of about 12 months.
}
\format{
\code{fisherB1} is a vector of 254 observations (in the format hours.minutes).
\code{fisherB1c} contains the same observations in a circular objects (minutes are expressed as decimals).
}

\source{
Cox, D.R. and Lewis, P.A.W. (1966) The Statistical Analysis of Series of Events. London : Methuen & CO. Ltd. pp. 254-255}

\seealso{
N.I. Fisher (1993) Statistical analysis of circular data. Cambridge University Press. Pag. 239.
}

\examples{
data(fisherB1c)
par(mfcol=c(1,2))
plot(fisherB1c, main="Clock 24", shrink=1.5)
plot(fisherB1c, template="clock12", main="Clock 12", shrink=1.5)
}

\keyword{datasets}
