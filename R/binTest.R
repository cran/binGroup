"binTest" <-
function(n, Y, p.hyp, alternative="two.sided", method="Exact")
{

 esti=Y/n 

 if(method=="Exact"){p.val=binom.test(x=Y,n=n,p=p.hyp, alternative=alternative)$p.value}
 
 if(method=="Score")
  {
   teststat = (esti-p.hyp)/(sqrt(p.hyp*(1-p.hyp)/n))
   if(alternative=="less"){p.val = pnorm(q=teststat,lower.tail=TRUE)}
   if(alternative=="greater"){p.val = pnorm(q=teststat,lower.tail=FALSE)} 
   if(alternative=="two.sided"){p.val= min( 2*pnorm(q=teststat, lower.tail = FALSE) , 2*pnorm(q=teststat, lower.tail = TRUE), 1)} 
  }

 if(method=="Wald")
  {
   teststat = (esti-p.hyp)/(sqrt(esti*(1-esti)/n))
   if(alternative=="less"){p.val = pnorm(q=teststat,lower.tail=TRUE)}
   if(alternative=="greater"){p.val = pnorm(q=teststat,lower.tail=FALSE)} 
   if(alternative=="two.sided"){p.val= min( 2*pnorm(q=teststat, lower.tail = FALSE) , 2*pnorm(q=teststat, lower.tail = TRUE), 1)}  
  }
 
out<-list(p.val=p.val,
estimate=esti,
alternative=alternative,
p.hyp=p.hyp)

class(out)<-"bgtTest"
out
}

