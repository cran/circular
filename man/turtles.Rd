\name{turtles}
\alias{turtles}
\docType{data}
\title{
Arrival directions of displaced sea turtles}
\description{
The _turtles_ dataset has 10 rows and 2 columns. The observations are the directions from which 10 green sea turtles approached their nesting island (Ascension Island, South Atlantic Ocean) after having been displaced to open-sea sites.  
}
\usage{data(turtles)}
\format{
  A data frame with 10 observations on the following 2 variables.
  \describe{
    \item{\code{id}}{a numeric vector: the turtle ID}
    \item{\code{arrival}}{a numeric vector: the direction of arrival to Ascension Island}
  }
}
\source{
Luschi, P., Akesson, S., Broderick, A. C., Glen, F., Godley, B. J., Papi F., and Hays, G. C. (2001) Testing the navigational abilities of ocean migrants: displacement experiments on green sea turtles (_Chelonia mydas_). Behav. Ecol. Sociobiol. (50):528-534.  
  
}
\examples{
data(turtles)
turtles[,2] <- circular(turtles[,2], units='degrees', template='geographics')
plot(turtles[,2])
}
\keyword{datasets}
