\name{fitmssg}
\alias{fitmssg}
\title{Computing the maximum likelihood estimator for the mixtures of skewed sub-Gaussian stable distributions using the EM algorithm.}
\description{Each \eqn{d}-dimensional skewed sub-Gaussian stable (SSG) random vector \eqn{\bf{Y}}, admits the representation given by Teimouri (2022):
\deqn{
{\bf{Y}} \mathop=\limits^d {\boldsymbol{\mu}}+\sqrt{P}{\boldsymbol{\lambda}}\vert{Z}_0\vert + \sqrt{P}{\Sigma}^{\frac{1}{2}}{\bf{Z}}_1,
}
where \eqn{\boldsymbol{\mu}} (location vector in \eqn{{{{R}}}^{d}}, \eqn{\boldsymbol{\lambda}} (skewness vector in \eqn{{{{R}}}^{d}}), \eqn{\Sigma} (positive definite symmetric dispersion matrix), and \eqn{0<\alpha \leq 2} (tail thickness) are model parameters. Furthermore, \eqn{P} is a positive stable random variable, \eqn{{Z}_0\sim N({0},1)}, and \eqn{\bf{Z}_1\sim N_{d}\bigl({\bf{0}}, \Sigma\bigr)}. We note that \eqn{Z}, \eqn{Z_0}, and \eqn{\boldsymbol{Z}_1} are mutually independent.
}

\usage{fitmssg(Y, K, eps = 0.15, initial = "FALSE", method = "moment", starts = starts)}
\arguments{
\item{Y}{an \eqn{n\times d} matrix of observations.}
\item{K}{number of component.}
\item{eps}{threshold value for stopping EM algorithm. It is 0.15 by default. The algorithm can be implemented faster if \code{eps} is larger.}
\item{initial}{logical statement. If \code{initial = TRUE}, then a list of the initial values must be given. Otherwise, it is determined by \code{method}.}
\item{method}{either \code{em} or \code{moment}. If \code{method = "moment"}, then the initial values are determined through the method of moment applied to each of \eqn{K} clusters that are obtained through the k-means method of Hartigan and Wong (1979). Otherwise, the initial values for each cluster are determined through the EM algorithm (Teimouri et al., 2018) developed for sub-Gaussian stable distributions applied to each of \eqn{K} clusters.}
\item{starts}{a list of initial values if \code{initial="TRUE"}. The list contains a vector of length \eqn{K} of mixing (weight) parameters, a vector of length \eqn{K} of tail thickness parameters, \eqn{K} vectors of length of \eqn{d} of location parameters, \eqn{K} dispersion matrices, \eqn{K} vectors of length of \eqn{d} of skewness parameters, respectively.}
}
%\details{}

\value{a list of estimated parameters corresponding to \eqn{K} clusters, predicted labels for clusters, the log-likelihood value across iterations, the Bayesian information criterion (BIC), and the Akaike information criterion (AIC).
}

\references{
M. Teimouri, 2022. Finite mixture of skewed sub-Gaussian stable distributions, arxiv.org/abs/2205.14067.

M. Teimouri, S. Rezakhah, and A. Mohammadpour, 2018. Parameter estimation using the EM algorithm for symmetric stable
  random variables and sub-Gaussian random vectors, \emph{Journal of Statistical Theory and Applications}, 17(3), 439-41.

J. A. Hartigan, M. A. Wong, 1979. Algorithm as 136: A k-means clustering algorithm, \emph{Journal of the Royal Statistical Society. Series c (Applied
  Statistics)}, 28, 100-108.
}
\author{Mahdi Teimouri}
\examples{
\donttest{
data(bankruptcy)
out1<-fitmssg(bankruptcy[,2:3], K=2, eps = 0.15, initial="FALSE", method="moment", starts=starts)
n1 <- 100
n2 <- 50
omega1 <- n1/(n1 + n2)
omega2 <- n2/(n1 + n2)
alpha1 <- 1.6
alpha2 <- 1.6
mu1 <- c(-1, -1)
mu2 <- c(6, 6)
sigma1 <- matrix( c(2, 0.20, 0.20, 0.5), 2, 2 )
sigma2 <- matrix( c(0.4, 0.10, 0.10, 0.2  ), 2, 2 )
lambda1 <- c(5, 5)
lambda2 <- c(-5, -5)
Sigma <- array( NA, c(2, 2, 2) )
Sigma[, , 1] <- sigma1
Sigma[, , 2] <- sigma2
starts<-list( c(omega1,omega2), c(alpha1,alpha2), rbind(mu1,mu2), Sigma, rbind(lambda1,lambda2) )
Y <- rbind( rssg(n1 , alpha1, mu1, sigma1, lambda1),  rssg(n2, alpha2, mu2, sigma2, lambda2) )
out2<-fitmssg(Y, K=2, eps=0.15, initial="TRUE", method="moment", starts=starts)
}
}
