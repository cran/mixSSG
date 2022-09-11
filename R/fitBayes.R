fitBayes <- function(y, mu0, sigma0, gamma0, delta0, epsilon)
{
	n <- length(y)
	N <- 3000
	M <- 1000
	n.i <- 5
	n.f <- 8
	n.burn <- 10
	tail.hat  <- rep(NA, n.f)
	N.slope <- 10
 	n.slope <- 5
	k <- 60
	K <- seq(1, k)
	z.y <- rz0 <- d0 <- NULL
	u.f <- rep(NA, 2)
	pdf <- rep(NA, n)
	theta<-matrix(NA, nrow = M, ncol = 3)  
	log.likelihood <- rep(NA, M)
	cri.seq <- c(1, seq(n.slope, M, n.slope) )
	f.alpha <- function(x, a, d) a*log(x) - x^a - d*x^2/2 + log( 1/2 )
    fprim.alpha <- function(x, a, d) a/x - a*x^(a - 1) - d*x 
	alpha.hat <- min(1.98,  1/sqrt( abs(6*var( log( abs( y[which( y != 0 )]) ) )/pi^2 - 1/2) ) )
	iq <- quantile(y, c(0.25, 0.75) )
	sigma.hat <- (iq[[2]] - iq[[1]])/ (1.46 - 0.0077*alpha.hat + 0.5501*alpha.hat^2)^(1/alpha.hat)
	mu.hat <- median(y)
	theta[1,] <- c(alpha.hat, sigma.hat, mu.hat)
	l <- 2
	cri <- 0.5
		while (cri < 1)
		{
		rz0 <- y - mu.hat
		d0 <- abs( rz0 )
		#vv <- alpha.hat + 2*sigma.hat*(exp(lgamma(alpha.hat*k/2+alpha.hat/2+1)+lgamma(alpha.hat*
		#k/2+alpha.hat/2+1/2)-lgamma(alpha.hat*k/2+1)-lgamma(alpha.hat*k/2+1/2))/((k+1)))^(1/alpha.hat)
			for (i in 1:n)
			{
 				#if (d0[i] > vv){
				#sum(1/(pi)*(-1)^(K-1)*exp(lgamma(a*K/2+1)-lgamma(K+1))*sin(K*pi*a*.5)*u^(-a*K/2-1))
				#pdf[i] <- 1/pi*sum( ( 2*sigma.hat/d0[i] )^(alpha.hat*K+1)*(-1)^(K-1)*exp(lgamma(
				#alpha.hat*K/2+1)+lgamma(alpha.hat*K/2+1/2)-lgamma(K+1))*sin(K*pi*alpha.hat*0.5) )
				#1/pi*sum( ( 2*sigma.hat/d0[i] )^(alpha.hat*K+1)*(-1)^(K-1)*exp(lgamma(alpha.hat*
				#K/2+1)+lgamma(alpha.hat*K/2+1/2)-lgamma(K+1))*sin(K*pi*alpha.hat*0.5) )
				#pdf[i] <- 1/pi*sum( ( 2*sigma.hat/d0[i] )^(alpha.hat*K-1)*(-1)^(K-1)*exp(lgamma(
				#alpha.hat*K/2+1)+lgamma(alpha.hat*K/2-1/2)-lgamma(K+1))*sin(K*pi*alpha.hat*0.5) )
				#z.y[i] <- t1/pdf[i]
				#z.y[i] <- 1/sqrt(4*sigma.hat^2*pi)*sum( (2*sigma.hat)^(alpha.hat*K+2)/(d0[i]^
				#(alpha.hat*K+3))*(-1)^(K-1)*exp(lgamma(alpha.hat*K/2+1)+
				#lgamma(alpha.hat*K/2+3/2)-lgamma(K+1))*sin(K*pi*alpha.hat*0.75) )/pdf[i]
 					#}else{
					#jj[i]<-i; if(!is.numeric(jj[i]))print(jj[i])
			r <- rpstable(N, alpha.hat); 
			pdf[i] <- 1/sqrt(4*sigma.hat^2*pi)*sum( exp(- 1/2*log(r) - d0[i]^2/(4*r*sigma.hat^2) ) )
			z.y[i] <- 1/sqrt(4*sigma.hat^2*pi)*sum( exp(- 3/2*log(r) - d0[i]^2/(4*r*sigma.hat^2) ) )/pdf[i]
						 #}
			}
		mu.hat <- ( sum(y*z.y)/2 + mu0*sigma.hat^2*sigma0^(-2) )/( sum(z.y)/2 + sigma.hat^2*sigma0^(-2) )
		sigma.hat <- sqrt(  ( sum( (y - mu.hat )^2*z.y )/2  + 2*delta0 )/( n + 2*gamma0 )   )
		obs.tr <- which( abs( rz0 ) <= 10e-18 )
		n.tr <- length( obs.tr )
		if( n.tr > 0 ) rz0 <- rz0[-obs.tr]
		n0 <- n - n.tr
			for (r1 in 1:n.f)
			{
				re0  <- sqrt( rexp(n0, 1) )
				Y.t0 <- rz0/re0
				dis0 <- (Y.t0/(sqrt(2)*sigma.hat))^2
				i0   <- 1
				m0 <- 0
				z0 <- rep(NA, n0)
				while(i0 <= n0)
				{
					f00 <- function(w) w^alpha.hat*exp( -w^alpha.hat - dis0[i0]*w^2/2 )/2
					if( dis0[i0] > 10e+20 )
					{
						z0[i0] <- 1/dis0[i0]
						}else{
							if( dis0[i0] < 0.1 & alpha.hat < 1 )
							{
								g   <- function(u) -alpha.hat + u^alpha.hat*alpha.hat + u^2*dis0[i0]
								w.0 <- uniroot( g, c(0, 10e+80) )$root
								g0  <- function(u) alpha.hat*log(u) -u^alpha.hat - u^2*dis0[i0]/2 + 1/2 + 10*log(10)
								w.end <- uniroot( g0, c(w.0, 10e+80) )$root
								w.old <- w.0
								j0 <- 1
									while(j0 <= 5)
									{
										f.u <- runif(1, 0, w.old^(alpha.hat)*exp( -w.old^alpha.hat - 
										dis0[i0]*w.old^2/2 )/2 )
										g00 <- function(w) w^(alpha.hat)*exp( -w^alpha.hat - dis0[i0]*w^2/2 )/2 - f.u
										u.f <- uniroot.all( g00 , c(0, w.end) )[1:2]
										if( f.u < 10e-16 ) u.f <- c(0, w.end)
										w.old <- ifelse( any( is.na(u.f) ) == TRUE,  w.0,  runif( 1, u.f[1], u.f[2] ) )
									j0 <- j0 + 1
									}
								z0[i0] <- w.old
								}else{
								z0[i0] <- ars(1, f.alpha, fprim.alpha, x = c( 0.01, 0.1, 1, 2, 5, 10, 40, 200, 2000,
										10e+5, 10e+8, 10e+30), m = 12, lb = TRUE, xlb = 0, a = alpha.hat, d = dis0[i0])
								}
						}
				i0 <- i0 + 1
				}
				f0 <- function(par) sum( -log( par[1] ) -par[1]*log(z0) + z0^par[1] )
				tail.hat[r1] <- suppressWarnings( optimize(f0, lower = 0.01, upper = 1.99)$minimum )
			}
				alpha.hat <- median( tail.hat[n.i:n.f] )
				theta[l, ] <- c(alpha.hat, sigma.hat, mu.hat)
				log.likelihood[l - 1] <- sum( log( sum(pdf) ) )
        		message("log.likelihood at iteration ", (l - 1), " is ", round( log.likelihood[l - 1], 5 ) )
				if( l >= N.slope & l %% n.slope == 0 )
				{
					i1   <- l/N.slope
					x1   <- log.likelihood[ (cri.seq[i1] + 1) : cri.seq[i1 + 1]       ]
					x2   <- log.likelihood[   cri.seq[i1 + 1] : (cri.seq[i1 + 2] - 1) ]
					x1.t <- subset( x1, x1 < quantile(x1, 0.8) & x1 > quantile(x1, 0.2) )
					x2.t <- subset( x2, x2 < quantile(x2, 0.8) & x2 > quantile(x2, 0.2) )
					slope1 <- lm( x1.t ~ seq(1:length(x1.t)) )$coefficients[[2]]
					slope2 <- lm( x2.t ~ seq(1:length(x2.t)) )$coefficients[[2]]
					if( abs( slope2 - slope1 ) <= epsilon && slope2 <= epsilon )	cri <- 1
				}
	l <- l + 1
	} 
			n.p <- 3
			aic <- -2*log.likelihood[l - 2] + 2*n.p
			bic <- -2*log.likelihood[l - 2] + n.p*log(n)
			theta.hat    <- colMeans( theta[(N.slope - 1):(l - 1), ] )
			names(alpha.hat) <- paste0("alpha")
			names(mu.hat)    <- paste0("mu")
			names(sigma.hat) <- paste0("sigma")
return( c( alpha.hat, sigma.hat, mu.hat, max.iter = l-1, list( log.likelihood = log.likelihood[ 1:(l-2) ],
 BIC = bic, AIC = aic  ) ) )
}
