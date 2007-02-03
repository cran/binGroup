"bgtSOC" <-
function(n,s,Y,conf.level=0.95,alternative="two.sided")

{
esti<-Y/n
kappa<-qnorm(conf.level)
eta<-(kappa^2)/3 + 1/6
gamma1<-((13/18)*kappa^2 + 17/18)*(-1)
gamma2<-(kappa^2)/18 + 7/36

midpo<-(Y+eta)/(n+2*eta)

if(alternative=="less")
  {upper = midpo + kappa * sqrt(esti*(1-esti) + (gamma1*esti*(1-esti) + gamma2)/n)/sqrt(n)
  CI=c( 0 ,upper)
  if(Y==n||upper>1){CI=c(0,1)}
  else{ CI=c( 0 ,upper)}
 }

if(alternative=="greater")
  {CI=c( midpo - kappa*sqrt(esti*(1-esti) + (gamma1*esti*(1-esti) + gamma2)/n)/sqrt(n) , 1)
  if(Y==0){CI=c(0,1)} }

if (alternative=="two.sided")
{
kappa<-qnorm(1-(1-conf.level)/2)
eta<-(kappa^2)/3 + 1/6
gamma1<-((13/18)*kappa^2 + 17/18)*(-1)
gamma2<-(kappa^2)/18 + 7/36

lower= midpo - kappa*sqrt(esti*(1-esti) + (gamma1*esti*(1-esti) + gamma2)/n)/sqrt(n)  
upper= midpo + kappa*sqrt(esti*(1-esti) + (gamma1*esti*(1-esti) + gamma2)/n)/sqrt(n)
 
if(Y==0){CI=c(0,upper)} 
else{if(Y==n||upper>1){CI=c(lower,1)}
else{CI=c(lower, upper)}}
}

KI=c( 1-(1-CI[1])^(1/s) ,  1-(1-CI[2])^(1/s))
KI

}

