dssg  <- function(Y, alpha, Mu, Sigma, Lambda)
{  # computes density function of the skewed sub-Gaussian distribution
	n     <- ifelse( is.null( dim(Y) ), 1 , length(Y[,1]) )
	Dim   <- length(Mu)
	np    <- round( min(160, 160/alpha) )
	n.p   <- seq(1, np)
	n.sim <- 5000
	pdf   <- rep(NA, n)
	dis   <- rep(NA, n)
	m     <- rep(NA, n)
     	lb    <- 2 + 2*(exp(lgamma(alpha*np/2+1+alpha/2) - lgamma(alpha*np/2+1) + lgamma((Dim+np*alpha)/2+alpha/2) -
			lgamma((Dim+np*alpha)/2))/(np+1))^(2/alpha)
	Omega <- Sigma + Lambda%*%t(Lambda)
	delta <- 1 - mahalanobis(Lambda, rep(0, Dim), Omega)
	C0    <- suppressWarnings( 2/( (2*pi)^(Dim/2+1/2)*sqrt(det(Sigma)) ) )
		for (i in 1:n)
		{
            	dis[i] <- ifelse( n == 1, mahalanobis(Y, Mu, Omega), mahalanobis(Y[i, ], Mu, Omega) )
			m[i] <-  ifelse( n == 1, (t(c(Lambda))%*%solve(Omega)%*%c(Y - Mu))[1],
					(t(c(Lambda))%*%solve(Omega)%*%c(Y[i, ] - Mu))[1] )
			if ( dis[i] >= lb )
			{
				pdf[i] <- sqrt(dis[i]*delta/pi)*sum( (-1)^(n.p-1)*sin(n.p*pi*alpha/2)*
						exp(lgamma(alpha*n.p/2+1) - lgamma(n.p+1) + lgamma((Dim+n.p*alpha)/2) )*
						(dis[i]/2)^(-(Dim+1+n.p*alpha)/2)*pt(q = m[i]*sqrt((Dim+n.p*alpha)/
						(dis[i]*delta)), df = Dim+n.p*alpha)  )
			}
			else
			{
              		r.p    <- rpstable(n.sim, alpha)
               		pdf[i] <- sqrt(2*pi*delta)*mean( exp(-Dim/2*log(r.p) - dis[i]/(2*r.p) + pnorm(m[i], 0, sqrt(r.p*delta), log.p=TRUE) ) )
            	}
		}
return( pdf )
}
