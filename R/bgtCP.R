"bgtCP" <-
function(n, Y, s, conf.level=0.95, alternative="two.sided")
{
lower<-0
upper<-1
if(alternative=="two.sided")
{
 if(Y!=0)
   {lower<-qbeta((1-conf.level)/2, Y, n-Y+1)}

 if(Y!=n)
   {upper<-qbeta(1-(1-conf.level)/2, Y+1, n-Y)}
}

if(alternative=="less")
{
 if(Y!=n)
   {upper<-qbeta(1-(1-conf.level), Y+1, n-Y)}
}

if(alternative=="greater")
{
 if(Y!=0)
   {lower<-qbeta((1-conf.level), Y, n-Y+1)}
}

estimate=1-(1-Y/n)^(1/s)


KI=c(1-(1-lower)^(1/s),1-(1-upper)^(1/s))

KI   
}

