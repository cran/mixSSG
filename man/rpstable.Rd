\name{rpstable}
\alias{rpstable}
\title{Simulating positive stable random variable.}
\description{The cumulative distribution function of positive stable distribution is given by
\deqn{
F_{P}(x)=\frac{1}{\pi}\int_{0}^{\pi}\exp\Bigl\{-x^{-\frac{\alpha}{2-\alpha}}a(\theta)\Bigr\}d\theta,
}
where \eqn{0<\alpha \leq 2} is tail thickness or index of stability and
\deqn{
a(\theta)=\frac{\sin\Bigl(\bigl(1-\frac{\alpha}{2}\bigr)\theta\Bigr)\Bigl[\sin \bigl(\frac{\alpha \theta}{2}\bigr)\Bigr]^{\frac{\alpha}{2-\alpha}}}{[\sin(\theta)]^{\frac{2}{2-\alpha}}}.
}
Kanter (1975) used the above integral transform to simulate positive stable random variable as
\deqn{P\mathop=\limits^d\Bigl(  \frac{a(\theta)}{W}  \Bigr)^{\frac{2-\alpha}{\alpha}},}
in which \eqn{\theta\sim U(0,\pi)} and \eqn{W} independently follows an exponential distribution with mean unity.
}
\usage{rpstable(n, alpha)}
\arguments{
\item{n}{the number of samples required.}
\item{alpha}{the tail thickness parameter.}
}
%\details{}

\value{
simulated realizations of size \eqn{n} from positive \eqn{\alpha}-stable distribution.
}
\references{
M. Kanter, 1975. Stable densities under change of scale and total variation inequalities, \emph{Annals of Probability}, 3(4), 697-707.
%J. M. Chambers, C. L. Mallows, and B. Stuck B, 1976. A method for simulating stable random variables. \emph{Journal of the American Statistical %Association}, 71, 354, 340-344.
}
\author{Mahdi Teimouri}
\examples{
\donttest{ rpstable(10, alpha = 1.2) }
}
