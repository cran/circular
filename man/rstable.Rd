\name{rstable}
\title{Random Generation from the Stable Family of Distributions}
\alias{rstable}

\description{
Returns random deviates from the stable family of probability distributions.
}

\usage{
rstable(n, scale = 1, index = stop("no index arg"), skewness = 0)
}

\arguments{
\item{n}{sample size.}
\item{index}{number from the interval (0, 2]. An index of 2 corresponds to the normal, 1 to the Cauchy.  Smaller values mean longer tails.}
\item{skewness}{number giving the modified skewness (see Chambers et al., 1976).  Negative values correspond to skewness to the left (the median is smaller than the mean, if it exists), and positive values correspond to skewness to the right (the median is larger than the mean). The absolute value of skewness should not exceed 1.}
\item{scale}{the scale of the distribution.}
}

\value{
random sample from the specified stable distribution.
}

\details{
This function return random variates from the Levy skew stable
distribution with \code{index}=$\alpha$, \code{scale}=$c$ and \code{skewness}=$\beta$.  The \code{skewness} parameter must lie in the range [-1,1] while the \code{index} parameter must lie in the range (0,2].  The Levy skew stable probability distribution is defined by a fourier transform,
\eqn{
p(x) = {1 \over 2 \pi} \int_{-\infty}^{+\infty} dt \exp(-it x - |c t|^\alpha (1-i \beta sign(t) \tan(\pi\alpha/2)))
}

When $\alpha = 1$ the term $\tan(\pi \alpha/2)$ is replaced by $-(2/\pi)\log|t|$.   For $\alpha = 2$ the distribution reduces to a Gaussian distribution with $\sigma = \sqrt{2} scale$ and the skewness parameter has no effect.  
For $\alpha < 1$ the tails of the distribution become extremely
wide.  The symmetric distribution corresponds to $\beta =0$.

The Levy alpha-stable distributions have the property that if $N$
alpha-stable variates are drawn from the distribution $p(c, \alpha,
\beta)$ then the sum $Y = X_1 + X_2 + \dots + X_N$ will also be
distributed as an alpha-stable variate, $p(N^{1/\alpha} c, \alpha, \beta)$.

There is no explicit solution for the form of $p(x)$ and there are no density, probability or quantile functions supplied for this distribution.
}

\author{Claudio Agostinelli}

\references{
Chambers, J. M., Mallows, C. L. and Stuck, B. W. (1976). A Method for Simulating Stable Random Variables. Journal of the American Statistical Association 71, 340-344.

Lo\"{g}ae\"{v}e, M. (1977). Probability Theory I. (fourth edition) Springer-Verlag, New York.
}

\seealso{
\link{rnorm}, \link{rcauchy}.
}

\examples{
hist(rstable(200, 1.5, .5)) #fairly long tails, skewed right
}

\keyword{distribution}
