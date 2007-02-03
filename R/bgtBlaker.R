"bgtBlaker" <-
function (n, Y, s, conf.level=0.95, alternative="two.sided")
{

# # # from the S code given in Blaker(2000), slightly changed # # #

tolerance=1e-04

acceptbin <- function(Y,n,p)
{
  p1 = 1-pbinom(Y-1, n, p)
  p2 = pbinom(Y, n, p)
  a1 = p1 + pbinom( qbinom(p1,n,p)-1, n, p )
  a2 = p2+1-pbinom( qbinom(1-p2,n,p), n, p )
  return(min(a1,a2))
}

lower<-0
upper<-1

if(Y!=0)
  {lower<-qbeta((1-conf.level)/2, Y, n-Y+1)
    {while(acceptbin(Y,n,lower+tolerance)<(1-conf.level))
    lower=lower+tolerance}
  }

if(Y!=n)
  {upper<-qbeta(1-(1-conf.level)/2, Y+1, n-Y)
    {while(acceptbin(Y,n,upper-tolerance)<(1-conf.level))
    upper=upper-tolerance}
  }
KI=c( 1-(1-lower)^(1/s), 1-(1-upper)^(1/s) )

KI  
}

