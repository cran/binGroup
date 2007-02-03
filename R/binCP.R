"binCP" <-
function(n, Y, conf.level=0.95, alternative="two.sided")
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
estimate=Y/n
conf.int=c(lower,upper)
conf.int
}

