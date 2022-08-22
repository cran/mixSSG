\name{rssg}
\alias{rssg}
\title{Simulating skewed sub-Gaussian \eqn{\alpha}-stable random vector.}
\description{Each skewed sub-Gaussian \eqn{\alpha}-stable (SSG) random vector \eqn{\bf{Y}}, admits the representation
\deqn{
{\bf{Y}} \mathop=\limits^d {\boldsymbol{\mu}}+\sqrt{P}{\boldsymbol{\lambda}}\vert{Z}_0\vert + \sqrt{P}{\Sigma}^{\frac{1}{2}}{\bf{Z}}_1,
}
where \eqn{{\boldsymbol{\mu}} \in {R}^{d} } is location vector, \eqn{{\boldsymbol{\lambda}} \in {R}^{d}} is skewness vector, \eqn{\Sigma} is a positive definite symmetric dispersion matrix, and \eqn{0<\alpha \leq 2} is tail thickness. Further, \eqn{P} is a positive \eqn{\alpha}-stable random variable, \eqn{{Z}_0\sim N({0},1)}, and \eqn{{\bf{Z}}_1\sim N_{d}\bigl({\bf{0}}, \Sigma\bigr)}. We note that \eqn{Z}, \eqn{Z_0}, and \eqn{{\bf{Z}}_1} are mutually independent.
}

\usage{rssg(n, alpha, Mu, Sigma, Lambda)}
\arguments{
\item{n}{the number of samples required.}
\item{alpha}{tail thickness parameter.}
\item{Mu}{a vector giving the location parameter.}
\item{Sigma}{a positive definite symmetric matrix specifying the dispersion matrix.}
\item{Lambda}{a vector giving the skewness parameter.}
}
%\details{}

\value{
simulated realizations of size \eqn{n} from the skewed sub-Gaussian \eqn{\alpha}-stable distribution.
}
%\references{}
\author{Mahdi Teimouri}
\examples{
n <- 4
alpha <- 1.4
Mu <- rep(0, 2)
Sigma <- diag(2)
Lambda <- rep(2, 2)
rssg(n, alpha, Mu, Sigma, Lambda)
}
