stoch <- function(Y, alpha0, Mu0, Sigma0, Lambda0)
{
  n.i <- 5
  n.f <- 8
  n.burn <- 10
  Dim <- length(Y[1, ])
  n <- length(Y[, 1])
  tail.hat  <- rep(NA, n.f)
  z0 <- NULL
  u.f <- rep(NA, 2)
  dis0 <- NULL
  m0 <- NULL
  tail.hat[1] <- alpha0
  Omega0 <- Sigma0 + Lambda0%*%t(Lambda0)
  delta0 <- abs( 1 - mahalanobis(Lambda0, rep(0, Dim), Omega0) )
  C0 <- 2*alpha0*sqrt( delta0/( (2*pi)^Dim*det(Sigma0) ) )
  f.alpha <- function(x, a, d, del, Dim, m0){(Dim + a - 1)*log(x) - x^a - d*x^2/2 +
      log( pnorm( abs(m0)*x/sqrt(del), 0, 1 ) )}
  fprim.alpha <- function(x, a, d, del, Dim, m0){ (Dim + a - 1)/x - a*x^(a - 1) -
      d*x + abs(m0)/sqrt(del)*dnorm( abs(m0)*
                                       x/sqrt(del) )/pnorm( abs(m0)*x/sqrt(del) ) }
  rz0 <- sweep( Y, 2, Mu0, "-")
  obs.tr <- which( rowSums( abs( rz0 ) ) <= 10e-18 )
  n.tr <- length( obs.tr )
  if( n.tr > 0 ) rz0 <- rz0[-obs.tr, ]
  n0 <- n - n.tr
  m0.p1 <- as.vector( t(c(Lambda0))%*%solve(Omega0) )
  for (r1 in 1:n.f)
  {
    re0  <- sqrt( rexp(n0, 1) )
    Y.t0 <- rz0/re0
    dis0 <- mahalanobis(Y.t0, rep(0, Dim), Omega0)
    i0   <- 1
    m0 <- as.vector( m0.p1%*%t(Y.t0) )
    while(i0 <= n0)
    {
      f00<-function(w)w^(Dim + alpha0 - 1)*exp( -w^alpha0 - dis0[i0]*w^2/2 )*
        pnorm( abs(m0[i0])/sqrt(delta0)*w )
      if( dis0[i0] > 10e+20 )
      {
        z0[i0] <- 1/dis0[i0]
      }else{
        if( dis0[i0] < 0.1 & alpha0 < 1 )
        {
          g   <- function(u) -Dim -alpha0 +1 +u^alpha0*alpha0 + u^2*dis0[i0]
          w.0 <- uniroot( g, c(0, 10e+80) )$root
          g0  <- function(u) (Dim +alpha0 -1)*log(u) -u^alpha0 - u^2*dis0[i0]/2 +
            pnorm( abs(m0[i0])/sqrt(delta0)*u )+ 10*log(10)
          w.end <- uniroot( g0, c(w.0, 10e+80) )$root
          w.old <- w.0
          j0 <- 1
          while(j0 <= 5)
          {
            f.u <- runif(1, 0, w.old^(Dim + alpha0 - 1)*exp( -w.old^alpha0 -
                                                               dis0[i0]*w.old^2/2 )*pnorm( abs(m0[i0])/sqrt(delta0)*w.old ) )
            g00 <- function(w){ w^( Dim + alpha0 - 1)*exp( -w^alpha0 -
                                                             dis0[i0]*w^2/2 )*pnorm( abs(m0[i0])/sqrt(delta0)*w ) - f.u }
            u.f <- uniroot.all( g00 , c(0, w.end) )[1:2]
            if( f.u < 10e-16 ) u.f <- c(0, w.end)
            w.old <- ifelse( any( is.na(u.f) ) == TRUE,  w.0,  runif( 1,
                                                                      u.f[1], u.f[2] ) )
            j0 <- j0 + 1
          }
          z0[i0] <- w.old
        }else{
          z0[i0] <- ars(1, f.alpha, fprim.alpha, x = c( 0.01, 0.1, 1, 2, 5, 10, 40,
                                                        200, 2000, 10e+5, 10e+8, 10e+30), m = 12, lb = TRUE, xlb = 0, a = alpha0,
                        d = dis0[i0], del = delta0, Dim = Dim, m0 = m0[i0])
        }
      }
      i0 <- i0 + 1
    }
    f0 <- function(par) sum( -log(par[1]) - (Dim + par[1] - 1)*log(z0) + z0^par[1] )
    alpha0 <- suppressWarnings( optimize(f0, lower = 0.01, upper = 1.99)$minimum )
    #print(alpha0)
    tail.hat[r1 + 1] <- alpha0
  }
  alpha0 <- median( tail.hat[n.i:n.f] )
  return( list(alphahat = alpha0) )
}
