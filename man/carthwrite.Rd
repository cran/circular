\name{Carthwrite}
\alias{dcarthwrite}
\title{Carthwrite's Power-of-Cosine Density Function}
\description{Density for the Carthwrite's power-of-cosine distribution.}
\usage{dcarthwrite(x, mu, psi)}

\arguments{
  \item{x}{a vector. The \code{x} and \code{q} objects are coerced to class \code{\link{circular}}.}
  \item{mu}{the location angular parameter. The object is coerced to class \code{\link{circular}}.}
  \item{psi}{the positive shape parameter.}
}

\details{The Carthwrite's power-of-cosine distribution has density
	\deqn{
		f(x)=\frac{2^{(1/\psi)-1} \Gamma^2((1/\psi)+1) (1+\cos(x-\mu))^{1/\psi}}
		{\pi\Gamma((2/\psi)+1)},
	}{%
		f(x)=[2^{(1/\psi)-1} \Gamma^2((1/\psi)+1) (1+\cos(x-\mu))^{1/\psi}] / [\pi \Gamma((2/\psi)+1)],
	}
	for \eqn{0 \le x < 2\pi}{0 <= x < 2\pi}.
}

\value{The density}
\references{Carthwrite, D.E. (1963). The use of directional spectra in studying the output of a wave recorder on a moving ship. Ocean Wave Spectra , 203-218.}

\author{Federico Rotolo}

\keyword{distribution}
\keyword{circle}
\keyword{circular}
\keyword{Carthwrite}
\keyword{cosine}
