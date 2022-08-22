#rpstable <- function(n, alpha)    # simulating positive stable r.v.
#{
#	theta <- runif(n, 0, pi)
#	t1 <-
#sin( (1-alpha/2)*theta )*( sin(alpha/2*theta) )^( alpha/(2-alpha) )/( sin(theta) )^( 1/(1-alpha/2) )
#	t2 <- -log( runif(n) )
#    (t1/t2)^( (2-alpha)/alpha )
#}
rpstable <- function(n, alpha)    # simulating positive stable r.v.
{
  w     <- -log(runif(n))
  theta <- pi * (runif(n) - 1/2)
  sigma <- (cos(pi*alpha/4))^(2/alpha)
  r     <- sin(alpha/2*(pi/2 + theta))/( cos(alpha*pi/4)*cos(theta) )^(2/alpha)*
    ( cos(alpha*pi/4 + (alpha/2 - 1)*theta)/w )^( (2 - alpha)/alpha )
  ( r - tan(pi*alpha/4) )*sigma + sigma*tan(pi*alpha/4)
}
