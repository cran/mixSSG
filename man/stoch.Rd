\name{stoch}
\alias{stoch}
\title{Estimating the tail index of the skewed sub-Gaussian stable distribution using the stochastic EM algorithm given that other parameters are known.}
\description{
Suppose \eqn{{\boldsymbol{Y}}_1,{\boldsymbol{Y}}_2, \cdots,{\boldsymbol{Y}}_n} are realizations following \eqn{d}-dimensional skewed sub-Gaussian stable distribution. Herein, we estimate the tail thickness parameter \eqn{0<\alpha \leq 2} when \eqn{\boldsymbol{\mu}} (location vector in \eqn{{{{R}}}^{d}}, \eqn{\boldsymbol{\lambda}} (skewness vector in \eqn{{{{R}}}^{d}}), and \eqn{\Sigma} (positive definite symmetric dispersion matrix are assumed to be known.
}
\usage{stoch(Y, alpha0, Mu0, Sigma0, Lambda0)}
\arguments{
\item{Y}{a vector (or an \eqn{n\times d} matrix) at which the density function is approximated.}
\item{alpha0}{initial value for the tail thickness parameter.}
\item{Mu0}{a vector giving the initial value for the location parameter.}
\item{Sigma0}{a positive definite symmetric matrix specifying the initial value for the dispersion matrix.}
\item{Lambda0}{a vector giving the initial value for the skewness parameter.}
}
\details{Here, we assume that parameters \eqn{{\boldsymbol{\mu}}}, \eqn{{\boldsymbol{\lambda}}}, and \eqn{\Sigma} are known and only the tail thickness parameter needs to be estimated.}
\value{ Estimated tail thickness parameter \eqn{\alpha}, of the skewed sub-Gaussian stable distribution.}
%\references{}
\author{Mahdi Teimouri}
\examples{
n <- 100
alpha <- 1.4
Mu <- rep(0, 2)
Sigma <- diag(2)
Lambda <- rep(2, 2)
Y <- rssg(n, alpha, Mu, Sigma, Lambda)
\donttest{stoch(Y, alpha, Mu, Sigma, Lambda)}
}
