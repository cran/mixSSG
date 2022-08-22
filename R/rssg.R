rssg<-function(n, alpha, Mu, Sigma, Lambda)
{
  Dim <- length(Mu)
	Y <- matrix(NA, nrow = n, ncol = Dim)
		for (i in 1:n)
		{
			Z <- rpstable(1, alpha)
			X <- mvrnorm( 1, mu = rep(0, Dim), Sigma = Sigma )
			u <- abs( rnorm(1) )
			Y[i, ] <- Mu + Lambda*u*sqrt(Z) + sqrt(Z)*X
		}
 Y
}
