\name{JonesPewsey}
\alias{djonespewsey}
\alias{jonespewsey}

\title{Jones and Pewsey Density Function}
\description{
  Density for the Jones and Pewsey circular distribution. 
}

\usage{djonespewsey(x, mu, kappa, psi)}

\arguments{
  \item{x}{a vector. The object is coerced to class \code{\link{circular}}.}
  \item{mu}{direction parameter of the distribution. The object is coerced to class \code{\link{circular}}.}
  \item{kappa}{non-negative concentration parameter of the distribution.}
  \item{psi}{real shape parameter.}
}

\details{The JonesPewsey distribution has density
	\deqn{
		f(x)=\frac{(\cosh(\kappa\psi) + \sinh(\kappa\psi)\cos(x-\mu))^{1/\psi}}
		{2\pi P_{1/\psi}(\cosh(\kappa\psi))},
	}{%
		f(x)=[(\cosh(\kappa\psi) + \sinh(\kappa\psi)\cos(x-\mu))^{1/\psi}] / [2\pi P_{1/\psi}(\cosh(\kappa\psi))],	
	}
	for \eqn{0 \le x < 2\pi}{0 <= x < 2 \pi}, where \eqn{P_{1/\psi}(\cdot)}{P_{1/\psi}(.)} is the associated Legendre function of the first kind, 		degree \eqn{1/\psi} and order 0.
}

\value{The density}

\references{Jones , M.C. and Pewsey, A. (2005). A family of simmetric distributions on the circle. J. Am. Statist. Assoc. 100, 1422-1428}
\author{Federico Rotolo}

\examples{
ff <- function(x) djonespewsey(x, mu=circular(4), kappa=1.8, psi=-.6)
curve.circular(ff, join=TRUE, xlim=c(-1, 1), ylim=c(-1.2, 1.2),
  main="Density of a JonesPewsey Distribution",
  xlab=expression(paste(mu,"=1.3",pi,", ",kappa,"=1.8, ",psi,"=-0.6"))
  )
}

\keyword{distribution}
\keyword{circle}
\keyword{circular}
