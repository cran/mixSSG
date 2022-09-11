\name{fitBayes}
\alias{fitBayes}
\title{Estimating parameters of the symmetric \eqn{\alpha}-stable (S\eqn{\alpha}S) distribution using Bayesian paradigm.}
\description{
Let \eqn{{{y}}_1,{{y}}_2, \cdots,{{y}}_n} are \eqn{n} realizations form S\eqn{\alpha}S distribution with parameters \eqn{\alpha, \sigma}, and \eqn{\mu}. Herein, we estimate parameters of symmetric univariate stable distribution within a Bayesian framework. We consider a uniform distribution for prior of tail thickness, that is \eqn{\alpha \sim U(0,2)}. The normal and inverse gamma conjugate priors are designated for \eqn{\mu} and \eqn{\sigma^2} with density functions given, respectively, by
\deqn{
\pi(\mu)=\frac{1}{\sqrt{2\pi}\sigma_{0}}\exp\Bigl\{-\frac{1}{2}\Bigl(\frac{\mu-\mu_0}{\sigma_0}\Bigr)^{2}\Bigr\},
}
and
\deqn{
\pi(\delta)= \delta_{0}^{\gamma_{0}}\delta^{-\gamma_0-1}\exp\Bigl\{-\frac{\delta_0}{\delta}\Bigr\},
}
where \eqn{\mu_0 \in R}, \eqn{\sigma_0>0}, \eqn{\delta=\sigma^2}, \eqn{\delta_0>0}, and \eqn{\gamma_0>0}.
}
\usage{fitBayes(y, mu0, sigma0, gamma0, delta0, epsilon)}
\arguments{
\item{y}{ vector of realizations that following S\eqn{\alpha}S distribution.}
\item{mu0}{the location hyperparameter corresponding to \eqn{\pi(\mu)}.}
\item{sigma0}{the standard deviation hyperparameter corresponding to \eqn{\pi(\mu)}.}
\item{gamma0}{the shape hyperparameter corresponding to \eqn{\pi(\delta)}.}
\item{delta0}{the rate hyperparameter corresponding to \eqn{\pi(\delta)}.}
\item{epsilon}{a positive small constant playing the role of threshold for stopping sampler.}
}
\value{ Estimated tail thickness, location, and scale parameters, number of iterations to attain convergence, the log-likelihood value across iterations, the Bayesian information criterion (BIC), and the Akaike information criterion (AIC).}
%\references{}
\author{Mahdi Teimouri}
\examples{
n <- 100
alpha <- 1.4
mu <- 0
sigma <- 1
y <- rnorm(n)
\donttest{fitBayes(y, mu0 = 0, sigma0 = 0.2, gamma0 = 10e-5, delta0 = 10e-5, epsilon = 0.005)}
}
