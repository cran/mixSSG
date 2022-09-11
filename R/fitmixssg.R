fitmssg <- function(Y, K, eps = 0.15, initial = "FALSE", method="moment", starts = starts)
{
	if( !is.matrix(Y) ) Y <- as.matrix(Y)
	if( any( is.nan(Y) ) || any( Y == Inf ) || any( is.na(Y) ) ) stop("Y must be a numeric matrix.", call. = FALSE)
	n    <- length(Y[, 1])
	Dim  <- length(Y[1, ])
	if(n <= 10*K) stop("The sample size is small.", call. = FALSE)
	if( K != round(K) || K <= 1) stop("The number of components must be an integer greater than one.", call. = FALSE)
	M <- 1000
	M0 <- 5
	N0 <- 3
	n.burn <- 20
	N.slope <- 10
	n.slope <- 5
	cri.seq <- c(1, seq(n.slope, M, n.slope) )
	n.sim   <- 3000
	alpha.stoch <- rep(NA, M0)
	z <- rep(NA, n)
	Mu.hat    <- array(0, c(Dim, K, M) )
	Sigma.hat <- array(0, c(Dim, (Dim*K), M) )
	alpha.hat <- matrix(0, ncol = K, nrow = M)
	omega.hat <- matrix(0, ncol = K, nrow = M)
	tau.hat <- E0 <-E1 <- E2 <- E3  <- matrix(0, ncol = K, nrow = n)
	Lambda.hat <- array(0, c(Dim, K, M) )
	log.likelihood <- rep(NA, M)
	omega <- alpha <- Mu <- Sigma <- Lambda <- vector("list", K)
		if(initial == "FALSE")
		{
			Y0    <- Y[ mahalanobis( Y, apply(Y, 2, mean), cov(Y))/(n-1) + 1/n < 2*Dim/n, ]
			clust <- kmeans(Y0, K)
			n.x   <- rep(NA, K)
			if(method == "moment")
			{
				for(k in 1:K)
				{
					omega.hat[1, k] <- sum(clust$cluster == k)/n
					x      <- Y[clust$cluster == k, ]
					n.x[k] <- ifelse( is.null( dim(x)  ), 1 , length(x[,1]) )
					if(n.x[k] == 1) stop("\r ----- You need to enter another set of the initial values -----", call. = FALSE)
					Mu.hat[, k, 1] <- apply(x, 2, median)
					x.z <- suppressWarnings( abs( sweep(x, 1, Mu.hat[, k, 1], "-") ) )
					a.hat <-0
					for(i in 1:Dim)a.hat <- a.hat +  1/sqrt( abs(6*var( log( x.z[which( x.z[, i] != 0 ), i] ))/pi^2 - 1/2) )
					a.hat  <- min( 1.99, a.hat/Dim )
					Sigma.initial <- matrix(0, nrow = Dim, ncol = Dim)
					Lambda.hat[, k, 1] <-  sapply( 1:Dim, function(i) sign( mean( (x[, i] - mean(x[, i]))^3 ) ) )
					Sigma0 <- diag( sapply( 1:Dim, function(i) 2*(-log( abs( mean(cos(x[, i])) )))^(2/a.hat) ) )
					Sigma.initial <- Sigma0
						for(i in 1:(Dim-1))
						{
							for(j in (i+1):Dim)
							{
								x.m.25  <- quantile(x[, i] - x[, j], 0.25)[[1]]
								x.m.75  <- quantile(x[, i] - x[, j], 0.75)[[1]]
								x.p.25  <- quantile(x[, i] + x[, j], 0.25)[[1]]
								x.p.75  <- quantile(x[, i] + x[, j], 0.75)[[1]]
								p1 <- (x.p.75 - x.p.25)/(1.46 - 0.0077*a.hat + 0.5501*a.hat^2)^(1/a.hat)
								p2 <- (x.m.75 - x.m.25)/(1.46 - 0.0077*a.hat + 0.5501*a.hat^2)^(1/a.hat)
								Sigma.initial[ i, j] <-  (p1^2 - p2^2)/2
								Sigma.initial[j, i] <- Sigma.initial[i, j]
							}
						}
					Sigma.hat[ , ( (k-1)*Dim + 1 ):(k*Dim), 1] <- Sigma.initial
					if( any( eigen(Sigma.initial)$value < 0 )) Sigma.hat[ , ( (k-1)*Dim + 1 ):(k*Dim), 1] <- Sigma0
					alpha.hat[1, k] <- a.hat
				}
					if(min(n.x) <= 2) stop("\r --- You need to enter another set of the initial values ---", call. = FALSE)
			}
			if(method == "em")
			{
				M0  <- 20
				n.k <- 50
				for(k in 1:K)
				{
					omega.hat[1, k] <- sum(clust$cluster == k)/n
					x      <- Y[clust$cluster == k, ]
					n.x[k] <- ifelse( is.null( dim(x)  ), 1 , length(x[,1]) )
					if(n.x[k] == 1) stop("\r ----- You need to enter another set of the initial values -----", call. = FALSE)
					x.z <- suppressWarnings( abs( sweep(x, 1, apply(x, 2, median), "-") ) )
					Lambda.hat[, k, 1] <-  sapply( 1:Dim, function(i) sign( mean( (x[, i] - mean(x[, i]))^3 ) ) )
					a.hat <- 0
					for(i in 1:Dim) a.hat <- a.hat +  1/sqrt( abs(6*var( log( x.z[which( x.z[, i] != 0 ), i] ))/pi^2 - 1/2) )
					a.hat  <- min( 1.99, a.hat/Dim )
					Sigma0 <- cov(x)
					Mu0  <- apply(x, 2, median)
					e1 <- dis <- rep(NA, n.x[k])
					for (r in 1:M0)
					{
						n.s <- seq( 1, round( min(80, 80/a.hat) ) )
						L <- 2 + (exp(lgamma(a.hat*n.k/2 + a.hat/2 +1 ) + lgamma(a.hat*n.k/2 + a.hat/2 + Dim/2 + 1)-
								lgamma(a.hat*n.k/2 + 1) - lgamma(a.hat*n.k/2 + Dim/2 + 1))/( (n.k + 1) ) )^(2/a.hat)
						dis <- mahalanobis(x, Mu0, Sigma0)
							for (i in 1:n.x[k])
							{
								if (dis[i] > L)
								{
									s1 <- sum((-1)^n.s*dis[i]^(-a.hat*n.s/2-Dim/2-1)*
									exp(lgamma(a.hat*n.s/2+1)+lgamma(a.hat*n.s/2+Dim/2+1)-
									lgamma(n.s+1))*sin(n.s*pi*a.hat*.75) )
									s2 <- sum((-1)^n.s*dis[i]^(-a.hat*n.s/2-Dim/2)*
									exp(lgamma(a.hat*n.s/2+1)+lgamma(a.hat*n.s/2+Dim/2)-
									lgamma(n.s+1))*sin(n.s*pi*a.hat*.75) )
									e1[i] <- s1/s2
								}
								else
								{
									r.p <- rpstable(n.sim, a.hat)
									e1[i] <- sum( na.omit(r.p^(-Dim/2-1)*exp(-0.5*dis[i]/r.p)) )/
										   sum( na.omit(r.p^(-Dim/2)*exp(-0.5*dis[i]/r.p)) )
								}
							}
						Mu0 <- colSums(e1*x)/sum(e1)
						Sigma0 <- matrix( rowSums( sapply(1:n.x[k], function(i) e1[i]*c(x[i, ] - Mu0)%o%c(x[i, ] - Mu0) )), nrow = Dim, ncol = Dim
)/n.x[k]
					}
					Mu.hat[ ,k, 1] <- Mu0
					Sigma.hat[ , ( (k-1)*Dim + 1 ):(k*Dim), 1] <- Sigma0
					alpha.hat[1, k] <- a.hat
				}
					if(min(n.x) <= 2) stop("\r --- You need to enter another set of the initial values ---", call. = FALSE)
			}
		}
		if(initial == "TRUE")
		{
			if( any( starts[[1]] <=0 ) ) stop("The mixing proportions must be positive.", call. = FALSE)
			if( any( is.nan( starts[[1]] ) ) || any( starts[[1]] == Inf) || any( is.na( starts[[1]] ) ) || sum( starts[[1]] ) != 1) stop("Sum of
mixing proportions must be one.", call. = FALSE)
			if( any( is.nan( starts[[2]] ) ) || any( starts[[2]] == Inf) || any( is.na( starts[[2]] ) ) || any( starts[[2]] <= 0 ) || any(
starts[[2]] > 2 ) ) stop("The tail index must be in (0, 2].", call. = FALSE)
				#for(k in 1:K)
				#{
				#	omega.hat[1, k] <- starts[[1]][k]
				#	alpha.hat[1, k] <- starts[[2]][k]
				#	Mu.hat[, k, 1]  <- starts[[3 + 3*(k - 1) ]]
				#	Sigma.hat[, ( (k-1)*Dim + 1 ):(k*Dim), 1] <- starts[[4 + 3*(k - 1) ]]
				#	Lambda.hat[, k , 1] <- starts[[5 + 3*(k - 1)]]
				#}
		    for(k in 1:K)
		    {
		      omega.hat[1, k]    <- starts[[1]][[k]]
		      alpha.hat[1, k]    <- starts[[2]][[k]]
		      Mu.hat[, k, 1]     <- starts[[3]][k, ]
		      Sigma.hat[, , 1]   <- starts[[4]][, , k]
		      Lambda.hat[, k, 1] <- starts[[5]][k, ]
		    }
		}
  r <- 2
cri <- 0.5
while (cri < 1)
{
	likelihood <- 0
	weighthed.pdf <- Sum0 <- Sum1 <- Sum2 <- Sum21 <- Sum3 <- matrix(NA, nrow = n, ncol = K);
	C0 <- delta <- rep(0,K)
	my <- dis <- matrix(0, nrow= n, ncol= K);
	ki <- seq(1,80)
	kk <- 80
	for (i in 1:n)
	{
#print(i)
		for (k in 1:K)
		{
			alpha <- alpha.hat[(r-1), k]
			Omega <- Sigma.hat[, ((k-1)*Dim + 1):(k*Dim), (r-1)] + Lambda.hat[, k, (r-1)]%*%t(Lambda.hat[, k, (r-1)])
			delta[k] <- 1 - mahalanobis(Lambda.hat[, k, (r-1)], rep(0, Dim), Omega)
			my[i,k]  <- as.numeric(t(c(Lambda.hat[, k, (r-1)]))%*%solve(Omega)%*%c(Y[i,] - Mu.hat[, k, (r-1)]))
			dis[i,k] <- mahalanobis(Y[i,], Mu.hat[, k, (r-1)], Omega)
			C0[k] <- suppressWarnings( 2/( (2*pi)^(Dim/2+1/2)*sqrt(det(Sigma.hat[ , ((k-1)*Dim + 1):(k*Dim), (r-1)])) ) )
			ii <- 0
			lb0 <- 2 + 2*(exp(lgamma(alpha*kk/2+1+alpha/2) - lgamma(alpha*kk/2+1) +
lgamma((Dim+kk*alpha+2*ii)/2+alpha/2)-lgamma((Dim+kk*alpha+2*ii)/2))/(kk+1))^(2/alpha)
			r.p <- rpstable(n.sim, alpha)
		if( dis[i,k] > lb0 )
		{
			Sum0[i, k] <- sqrt(dis[i,k]*delta[k]/pi)*sum( (-1)^(ki-1)*sin(ki*pi*alpha/2)*
			exp(lgamma(alpha*ki/2+1) - lgamma(ki+1) + lgamma((Dim+ki*alpha+2*ii)/2) )*
			(dis[i,k]/2)^(-(Dim+1+ki*alpha+2*ii)/2)*pt(q = my[i,k]*sqrt((Dim+ki*alpha+2*ii)/
			(dis[i,k]*delta[k])), df = Dim+ki*alpha+2*ii)  )
		}
		else
		{
Sum0[i, k] <- mean( sqrt(2*pi*delta[k])*exp(  -Dim/2*log(r.p)-dis[i,k]/(2*r.p) + pnorm(my[i,k], 0, sqrt(r.p*delta[k]), log.p=TRUE) ) )
		}
			ii <- 1
			lb1 <- 2 + 2*(exp(lgamma(alpha*kk/2+1+alpha/2) - lgamma(alpha*kk/2+1) +
lgamma((Dim+kk*alpha+2*ii)/2+alpha/2)-lgamma((Dim+kk*alpha+2*ii)/2))/(kk+1))^(2/alpha)
		if( dis[i, k] > lb1 )
		{
			Sum1[i, k] <- sqrt(dis[i,k]*delta[k]/pi)*sum( (-1)^(ki-1)*sin(ki*pi*alpha/2)*
			exp(lgamma(alpha*ki/2+1) - lgamma(ki+1) + lgamma((Dim+ki*alpha+2*ii)/2) )*
			(dis[i,k]/2)^(-(Dim+1+ki*alpha+2*ii)/2)*pt(q = my[i,k]*sqrt((Dim+ki*alpha+2*ii)/
			(dis[i,k]*delta[k])), df = Dim+ki*alpha+2*ii)  )
		}
		else
		{
Sum1[i, k] <- mean(sqrt(2*pi*delta[k])*exp(-(Dim+2)/2*log(r.p)-dis[i,k]/(2*r.p) + pnorm(my[i,k], 0, sqrt(r.p*delta[k] ), log.p=TRUE) ) )
		}
lb21 <- 2 + 2*(exp(lgamma(alpha*kk/2+1+alpha/2) - lgamma(alpha*kk/2+1) + lgamma((Dim+1+kk*alpha)/2+alpha/2) -
		lgamma((Dim+1+kk*alpha)/2))/(kk+1))^(2/alpha)
		if( dis[i, k] > lb21 )
		{
			Sum21[i, k] <- delta[k]/pi*sum( (-1)^(ki-1)*exp( lgamma(alpha*ki/2+1) - lgamma(ki+1) +
						lgamma((Dim+1+ki*alpha)/2) )*sin(ki*pi*alpha/2)*( dis[i,k]/2 +
						my[i,k]^2/(2*delta[k]) )^( -(Dim+1+ki*alpha)/2 ) )
			Sum2[i, k] <- Sum21[i, k] + my[i,k]*Sum1[i,k]
		}
		else
		{
Sum21[i, k] <- delta[k]*mean( exp(-(Dim+1)/2*log(r.p) - dis[i,k]/(2*r.p) - my[i,k]^2/(2*delta[k]*r.p)) )
Sum2[i, k]  <- Sum21[i, k] + my[i,k]*Sum1[i,k]
		}

			lb3 <- 2 + 2*(exp( lgamma(alpha*kk/2+1+alpha/2) - lgamma(alpha*kk/2+1) +
				lgamma((Dim+2+kk*alpha)/2+alpha/2) - lgamma((Dim+2+kk*alpha)/2) )/(kk+1))^(2/alpha)
		if( dis[i,k] > lb3 )
		{
			nu <- Dim + 2 + ki*alpha
			b <- my[i,k]*sqrt(nu/(delta[k]*dis[i,k]))
Sum3[i, k] <- (delta[k]*dis[i,k])^(3/2)/sqrt(pi)*sum( (-1)^(ki-1)*sin(ki*pi*alpha/2)*
				exp( lgamma(alpha*ki/2+1) + lgamma(nu/2) - lgamma(ki+1) )/nu*(dis[i,k]/2)^(-(nu+1)/2)*(
				nu*(nu-1)/(nu-2)*pt(b*sqrt((nu-2)/nu), df = (nu-2) )/pt(b, df = nu) - nu )*pt(b, df = nu) ) +
				2*my[i,k]*Sum21[i, k] + my[i,k]^2*Sum1[i, k]
		}
		else
		{
Sum3[i, k] <- gamma(3/2)*sqrt(2)*delta[k]^(3/2)*mean( exp(-Dim/2*log(r.p)-dis[i,k]/(2*r.p))*
			(1 + sign(my[i,k])*pgamma(my[i,k]^2/(2*r.p*delta[k]), 3/2, 1)) ) +
			2*my[i,k]*Sum21[i, k] + my[i,k]^2*Sum1[i, k]
		}
		weighthed.pdf[i, k] <- omega.hat[(r-1), k]*C0[k]*Sum0[i,k]
		}
		tau.hat[i, ] <- weighthed.pdf[i, ]/sum(weighthed.pdf[i, ])
	}
omega.hat[r, ] <- colMeans(tau.hat)

		E1 <- tau.hat*Sum1/Sum0
		E2 <- tau.hat*Sum2/Sum0
		E3 <- tau.hat*Sum3/Sum0
			for (k in 1:K)
			{
				S1 <- sapply(1:n, function(i)E1[i, k]*c(Y[i, ] - Mu.hat[, k, (r-1)])%o%c(Y[i, ] - Mu.hat[, k, (r-1)]) )
				S2 <- sapply(1:n, function(i)E3[i, k]*Lambda.hat[, k, (r-1)]%o%Lambda.hat[, k, (r-1)] )
				S3 <- sapply(1:n, function(i)E2[i, k]*Lambda.hat[, k, (r-1)]%o%t(c(Y[i, ] - Mu.hat[, k, (r-1)])) )
				S4 <- sapply(1:n, function(i)E2[i, k]*c(Y[i, ] - Mu.hat[, k, (r-1)])%o%t(Lambda.hat[, k, (r-1)]) )
				T1 <- matrix( rowSums(S1 + S2 - S3 - S4), nrow = Dim, ncol = Dim )
				Mu.hat[, k, r]  <- ( colSums( E1[, k]*Y) - sum(E2[, k])*Lambda.hat[, k, (r-1)] )/sum(E1[, k])
				Lambda.hat[, k, r]  <- colSums( E2[, k]*sweep(Y, 2, Mu.hat[, k, r], "-") )/sum(E3[, k])
				Sigma.hat[, ((k-1)*Dim + 1):(k*Dim), r] <- T1/sum( tau.hat[, k] )
				log.likelihood[r - 1] <- sum( log( rowSums(weighthed.pdf) ) )
				cluster <- apply(tau.hat, 1, which.max)
				if( length(table(cluster)) == 1 || table(cluster)[[k]] <= Dim) stop("\r ----- You need to enter another set of the initial values
-----", call. = FALSE)
				alpha.hat[r, k] <- stoch( Y[cluster == k,], alpha.hat[(r-1), k], Mu.hat[, k, r], Sigma.hat[, ((k-1)*Dim + 1):(k*Dim), r], Lambda.hat[, k, r] )$alphahat
				}
        	#	cat("log.likelihood at iteration", (r-1), "is", log.likelihood[r-1], "\n")
        		message("log.likelihood at iteration ", (r-1), " is ", round( log.likelihood[r-1], 5 ) )
				if( r >= N.slope & r %% n.slope == 0 )
				{
					i1   <- r/N.slope
					x1   <- log.likelihood[ (cri.seq[i1] + 1) : cri.seq[i1 + 1] ]
					x2   <- log.likelihood[   cri.seq[i1 + 1] : (cri.seq[i1 + 2] - 1) ]
					x1.t <- subset( x1, x1 < quantile(x1, 0.8) & x1 > quantile(x1, 0.2) )
					x2.t <- subset( x2, x2 < quantile(x2, 0.8) & x2 > quantile(x2, 0.2) )
					slope1 <- lm( x1.t ~ seq(1:length(x1.t)) )$coefficients[[2]]
					slope2 <- lm( x2.t ~ seq(1:length(x2.t)) )$coefficients[[2]]
					if( abs( slope2 - slope1 ) <= eps && slope2 <= eps )	cri <- 1
				}
	r <- r + 1
}
			n.p <- 2*K*Dim + (Dim - 1) + K*( Dim*(Dim - 1)/2 + Dim)
			aic <- -2*log.likelihood[r - 2] + 2*n.p
			bic <- -2*log.likelihood[r - 2] + n.p*log(n)
			for(k in 1:K)
			{
				omega[[k]] <- apply( omega.hat[(N.slope-1):(r-1), ], 2, mean )[k]
				alpha[[k]] <- apply( alpha.hat[(N.slope-1):(r-1), ], 2, mean )[k]
				Mu[[k]]    <- apply( Mu.hat[ , , (N.slope-1):(r-1)], c(1:2), mean )[ , k]
				Sigma[[k]] <- apply( Sigma.hat[ ,( (k-1)*Dim + 1 ):(k*Dim), (N.slope-1):(r-1)], c(1:2), mean )
				Lambda[[k]]<- apply( Lambda.hat[ , , (N.slope-1):(r-1)], c(1:2), mean )[ , k]
				names(omega)[k] <- paste0("omega", k)
				names(alpha)[k] <- paste0("alpha", k)
				names(Mu)[k] <- paste0("Mu", k)
				names(Sigma)[k] <- paste0("Sigma", k)
				names(Lambda)[k] <- paste0("Lambda", k)
			}
#print(Mu.hat[, , r])
return( c( omega, alpha, Mu, Sigma, Lambda, list( cluster = cluster, log.likelihood = log.likelihood[ 1:(r - 2) ], BIC = bic, AIC = aic  ) ) )
}
