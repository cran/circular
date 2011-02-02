\name{KatoJones}
\alias{rkatojones}
\alias{dkatojones}

\title{Kato and Jones Density Function}
\description{Density and random generation for the Kato and Jones distribution.}

\usage{
rkatojones(n, mu, nu, r, kappa, control.circular=list())
dkatojones(x, mu, nu, r, kappa)
}

\arguments{
  \item{x}{the angular value the density must be computed in.}
  \item{n}{number of observations.}  
  \item{mu}{the Mobius 'mu' parameter. The object is coerced to class \code{\link{circular}}.}
  \item{nu}{the Mobius 'nu' parameter. The object is coerced to class \code{\link{circular}}.}
  \item{r}{the Mobius 'r' parameter. It must be in [0,1).}
  \item{kappa}{the positive vonMises parameter.}
  \item{control.circular}{the attribute of the resulting object.}  
}

\details{The Kato and Jones distribution has density
	\deqn{
		f(x)= \frac{1-r^2}{2\pi\mathcal I_0(\kappa)}
		\exp\left[
			\frac{\kappa\{ \xi\cos(x-\eta)-2r\cos\nu \}}
			{1+r^2-2r\cos(x -\gamma)}
		\right]\\
		\phantom{\exp[]} \times \frac1{1+r^2-2r\cos(x -\gamma)},
	}{%
	f(x)= [1-r^2]/[{2\pi\mathcal I_0(\kappa)}{1+r^2-2r\cos(x -\gamma)}]
			exp[{\kappa\{ \xi\cos(x-\eta)-2r\cos\nu \}} / {1+r^2-2r\cos(x -\gamma)}],
	}
	for \eqn{0 \le x < 2\pi}{0 <= x < 2 \pi}, 
	where \eqn{\gamma=\mu+\nu}, \eqn{\xi=\{r^4+2r^2\cos(2\nu)+1\}^{1/2}}{\xi={r^4+2r^2 cos(2\nu)+1}^{1/2}}
	and \eqn{\eta=\mu+\arg[ r^2\{\cos(2\nu)+i\sin(2\nu)\}+1 ]}{\eta=\mu+arg[ r^2{cos(2\nu)+i sin(2\nu)}+1 ]}.
	
	Original code for random generation is by Kato, S. and Jones, M.C. and can be found at the address 
	http://pubs.amstat.org/doi/suppl/10.1198/jasa.2009.tm08313/suppl_file/t08-313code.txt.
}

\value{The density.
\code{dkatojones} gives the density and \code{rkatojones} generates random deviates.

}

\references{Kato , S. and Jones, M.C. (2010). A family of distributions on the circle with links to, and applications arising from,
   Mobius transformation. J. Am. Statist. Assoc. 105, 249-262.}
\author{Federico Rotolo}

\examples{
data1 <- rkatojones(n=100, mu=circular(0), nu=circular(pi/4), r=.2, kappa=1)
plot(data1)

data1 <- rkatojones(n=100, mu=circular(pi/3), nu=circular(pi), r=.7, kappa=2.3)
plot(data1)

ff <- function(x) dkatojones(x, mu=circular(pi/3), nu=circular(pi), r=.7, kappa=2.3)
curve.circular(ff, join=TRUE, xlim=c(-1, 1), ylim=c(-1.2, 1.2),
  main="Density of a KatoJones Distribution",
  xlab=expression(paste(mu,"=",pi,"/3, ",nu,"=",pi,", r=0.7, ",kappa,"=2.3"))
  )
}

\keyword{distribution}
\keyword{circle}
\keyword{circular}
