"bgtTest" <-
function(n, Y, s, p.hyp, alternative="two.sided", method="Exact")
{

if( n<1 || length(n)!=1){stop("number of groups n must be specified as a single integer>0")}
if( s<1 || length(s)!=1){stop("group size s must be specified as a single integer>0")}
if( Y<0 || length(Y)!=1){stop("observed number of positive groups Y must be specified as a single integer>0")}
if(Y>n) {stop("number of positive tests Y can not be greater than number of groups n")}
if(p.hyp<0 || p.hyp>1 || length(p.hyp)!= 1){stop("the proportion in the hypothesis p.hyp must be specified as a single value between 0 and 1")}
if(alternative!="less" && alternative!="greater"&& alternative!="two.sided"){stop("argument alternative mis-specified, must be either 'less', 'greater' or 'two.sided'")}
if(method!="Exact" && method!="Score"&& method!="Wald"){stop("argument alternative mis-specified")}

estimate=1-(1-Y/n)^(1/s)

bgtsumProb<-function(x, n, s, p)
{
sumprob=0
for(i in x)
  {
  sumprob = sumprob + choose(n=n,k=i) * ((1-(1-p)^s)^(i)) * (1-p)^(s*(n-i))
  }
sumprob
}

if (method=="Exact")
{
 if(alternative=="less")
  {p.val = bgtsumProb(x=0:Y, n=n, s=s, p=p.hyp)}

 if(alternative=="greater")
  {p.val = bgtsumProb(x=Y:n, n=n, s=s, p=p.hyp)}

 if (alternative=="two.sided")
  {p.val = min( 2*(bgtsumProb(x=0:Y, n=n, s=s, p=p.hyp)), 2*(p.val = bgtsumProb(x=Y:n, n=n, s=s, p=p.hyp)), 1)}

}


if (method=="Wald")
{
esti = 1-(1-Y/n)^(1/s)

varesti = (1-(1-esti)^s)/(n*(s^2)*(1-esti)^(s-2))  # variance estimator (see Swallow, 1985)

teststat = (esti-p.hyp)/sqrt(varesti) 


 if(alternative=="less"){p.val = pnorm(q=teststat,lower.tail=TRUE)}

 if(alternative=="greater"){p.val = pnorm(q=teststat,lower.tail=FALSE)} 

 if(alternative=="two.sided"){p.val= min( 2*pnorm(q=teststat, lower.tail = FALSE) , 2*pnorm(q=teststat, lower.tail = TRUE), 1)} 

}

if(method=="Score")
{
esti = Y/n
t.hyp = 1-(1-p.hyp)^s
teststat = (esti-t.hyp)/(sqrt(t.hyp*(1-t.hyp)/n))

 if(alternative=="less"){p.val = pnorm(q=teststat,lower.tail=TRUE)}

 if(alternative=="greater"){p.val = pnorm(q=teststat,lower.tail=FALSE)} 

 if(alternative=="two.sided"){p.val= min( 2*pnorm(q=teststat, lower.tail = FALSE) , 2*pnorm(q=teststat, lower.tail = TRUE), 1)} 
}

out<-list(p.val=p.val,
estimate=estimate,
alternative=alternative,
p.hyp=p.hyp)

class(out)<-"bgtTest"
out
}

