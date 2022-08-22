\name{dssg}
\alias{dssg}
\title{Approximating the density function of skewed sub-Gaussian \eqn{\alpha}-stable distribution.}
\description{Suppose \eqn{d}-dimensional random vector \eqn{\boldsymbol{Y}} follows a skewed sub-Gaussian \eqn{\alpha}-stable distribution with density function
\eqn{f_{\boldsymbol{Y}}(\boldsymbol{y} | \boldsymbol{\Theta})} for \eqn{{\boldsymbol{\Theta}}=(\alpha,\boldsymbol{\mu},\Sigma, \boldsymbol{\lambda})} where \eqn{\alpha}, \eqn{\boldsymbol{\mu}}, \eqn{\Sigma}, and \eqn{\boldsymbol{\lambda}} are tail thickness, location, dispersion matrix, and skewness parameters, respectively. Herein, we give a good approximation for \eqn{f_{\boldsymbol{Y}}(\boldsymbol{y} | \boldsymbol{\Theta})}. First , for \eqn{{\cal{N}}=50}, define
\deqn{
L=\frac{\Gamma(\frac{{\cal{N}}\alpha}{2}+1+\frac{\alpha}{2})\Gamma\bigl(\frac{d+{\cal{N}}\alpha}{2}+\frac{\alpha}{2}\bigr)}{
\Gamma(\frac{{\cal{N}}\alpha}{2}+1)\Gamma\bigl(\frac{d+{\cal{N}}\alpha}{2}\bigr)({\cal{N}}+1)}.
}
If \eqn{d(\boldsymbol{y})\leq 2L^{\frac{2}{\alpha}}}, then
\deqn{
f_{\boldsymbol{Y}}(\boldsymbol{y} | {\boldsymbol{\Theta}}) \simeq
\frac{{C}_{0}\sqrt{2\pi \delta }}{N} \sum_{i=1}^{N} \exp\Bigl\{-\frac{d(\boldsymbol{y})}{2p_{i}}\Bigr\}\Phi \bigl( m| 0, \sqrt{\delta p_{i}} \bigr)p_{i}^{-\frac{d}{2}},
}
where, \eqn{p_1,p_2,\cdots, p_N} (for \eqn{N=3000}) are independent realizations following positive \eqn{\alpha}-stable distribution that are generated using command \code{rpstable(3000, alpha)}. Otherwise, if \eqn{d(\boldsymbol{y})> 2L^{\frac{2}{\alpha}}}, we have
\deqn{
f_{\boldsymbol{Y}}(\boldsymbol{y} | \boldsymbol{\Theta})\simeq
\frac{{C}_{0}\sqrt{d(\boldsymbol{y})\delta}}{\sqrt{\pi}}
\sum_{j=1}^{{\cal{N}}}\frac{ (-1)^{j-1}\Gamma(\frac{j\alpha}{2}+1)\sin \bigl(\frac{j\pi \alpha}{2}\bigr)}
{\Gamma(j+1)\bigl[\frac{d(\boldsymbol{y})}{2}\bigr]^{\frac{d+1+j\alpha}{2}}}\Gamma\Bigl(\frac{d+j\alpha}{2}\Bigr)
 T_{d+j\alpha}\biggl(m\sqrt{\frac{d+j\alpha}{d(\boldsymbol{y})\delta}}\biggr),
}
where \eqn{T_{\nu}(x)} and \eqn{\Phi(x|a,b)} denote the distribution function of the Student's \eqn{t} (distribution )with \eqn{\nu} degrees of freedom) and normal (mean \eqn{a} and standard deviation \eqn{b}) distributions at point \eqn{x}, respectively, and
%{C}_{0}=2 (2\pi)^{-\frac{d+1}{2}}\vert{\Sigma}\vert^{-\frac{1}{2}}
\eqn{{C_{0}=2 (2\pi)^{-\frac{d+1}{2}}|{\Sigma}|^{-\frac{1}{2}},}}
\eqn{d(\boldsymbol{y})=(\boldsymbol{y}-\boldsymbol{\mu})^{'}{{\Omega}^{-1}}(\boldsymbol{y}-\boldsymbol{\mu}),}
\eqn{{m}=\boldsymbol{\lambda}^{'}{{\Omega}}^{-1}(\boldsymbol{y}-\boldsymbol{\mu}),}
\eqn{{\Omega}={\Sigma}+\boldsymbol{\lambda}\boldsymbol{\lambda}^{'},}
\eqn{{\delta}=1-\boldsymbol{\lambda}^{'}{\Omega}^{-1}\boldsymbol{\lambda}}.
}
\usage{dssg(Y, alpha, Mu, Sigma, Lambda)}
\arguments{
\item{Y}{a vector (or an \eqn{n\times d} matrix) at which the density function is approximated.}
\item{alpha}{tail thickness parameter.}
\item{Mu}{a vector giving the location parameter.}
\item{Sigma}{a positive definite symmetric matrix specifying the dispersion matrix.}
\item{Lambda}{a vector giving the skewness parameter.}
}
%\details{}

\value{
simulated realizations of size \eqn{n} from positive \eqn{\alpha}-stable distribution.
}
%\references{
%M. Kanter, 1975. Stable densities under change of scale and total variation inequalities, \emph{Annals of Probability}, 3(4), 697-707.
%G. S. Mudholkar and A. D. Huston, 2001. The epsilon-skew-normal distribution for analyzing near-normal data, \emph{Journal of Statistical Planning and
%Inference}, 83, 291-309.
%
%J. M. Chambers, C. L. Mallows, and B. Stuck, 1976. A method for simulating stable random variables. \emph{Journal of the American Statistical Associati%on}, 71, 354, 340-344.
%}
\author{Mahdi Teimouri}
\examples{
n <- 4
alpha <- 1.4
Mu <- rep(0, 2)
Sigma <- diag(2)
Lambda <- rep(2, 2)
Y <- rssg(n, alpha, Mu, Sigma, Lambda)
dssg(Y, alpha, Mu, Sigma, Lambda)
}
