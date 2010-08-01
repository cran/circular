\name{swallows}
\alias{swallows}
\docType{data}
\title{
Orientation of juvenile barn swallows
}
\description{
The _swallows_ dataset has 114 rows and 2 columns. The observations are the headings of juvenile barn swallows (_Hirundo rustica_) tested in orientation cages (Emlen funnels) during autumn migration under simulated overcast conditions. 
}
\usage{data(swallows)}
\format{
  A data frame with 114 observations on the following 2 variables.
  \describe{
    \item{\code{treatment}}{a factor with levels \code{control} (control
      group: local magnetic field) and 
      \code{shifted} (shifted magnetic field, magnetic North =
      geographical West)}
    \item{\code{heading}}{a numeric vector:  modal heading of each bird}
  }
}
\source{
Giunchi, D., and Baldaccini N. E. (2004) Orientation of juvenile barn swallows (Hirundo rustica) tested in Emlen funnels during autumn migration. Behav. Ecol. Sociobiol. (56):124-131. 
}
\examples{
data(swallows)

swallows <- split(swallows$heading, swallows$treatment)
swallows <- lapply(swallows, function(x) circular(x, units='degrees', template='geographics'))

plot(swallows[[1]])
points(swallows[[2]], col=2)
legend('topright', legend=c('control', 'shifted'), pch=c(19,19), col=c(1,2))
}
\keyword{datasets}
