\name{wind}
\alias{wind}
\docType{data}
\title{
Col De La Roa wind direction  
}
\description{
  In a place named "Col de la Roa" in the Italian Alps there is a
  meteorological station that records via data-logger several
  parameters. Measures are made every 15 minutes, in this dataset
  we report the wind direction recorded every day from January 29,
  2001 to March 31, 2001 from 3.00am to 4.00am included. Which
  means 5 observations every day for a total of 310 measures.    
}
\usage{data(wind)}
\format{
  This data frame contains one variables (wind direction) in radians.
}

\source{
  http://www.tesaf.unipd.it/SanVito/dati.htm
}
\references{
  C. Agostinelli (2007) Robust estimation for circular data, Computational
  Statistics and Data Analysis, 51(12), 5867-5875,
  doi = doi:10.1016/j.csda.2006.11.002
}
\examples{
  data(wind)
  wind <- circular(wind, template='geographics')
  par(mfcol=c(1,2))
  plot(wind)
  plot(density(wind, bw=40), main='')
}
\keyword{datasets}
