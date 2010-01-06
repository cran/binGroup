"bgtCP" <-
function(n, y, s, conf.level=0.95, alternative="two.sided")
{
lower<-0
upper<-1
if(alternative=="two.sided")
{
 if(y!=0)
   {lower<-qbeta((1-conf.level)/2, y, n-y+1)}

 if(y!=n)
   {upper<-qbeta(1-(1-conf.level)/2, y+1, n-y)}
}

if(alternative=="less")
{
 if(y!=n)
   {upper<-qbeta(1-(1-conf.level), y+1, n-y)}
}

if(alternative=="greater")
{
 if(y!=0)
   {lower<-qbeta((1-conf.level), y, n-y+1)}
}

estimate=1-(1-y/n)^(1/s)


CI=c(1-(1-lower)^(1/s),1-(1-upper)^(1/s))

CI   
}

