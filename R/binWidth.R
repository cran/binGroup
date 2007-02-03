"binWidth" <-
function(n, p, conf.level=0.95, alternative="two.sided", method="CP")
{

# indicator function for the CI length at a special event
# in one sided case: length is defined as absolute difference between estimator and confidence bound

 L.Ind.bin<-function(Y, n, p, conf.level, alternative, method)

 {

  if(method=="Wald"){int=binWald(Y=Y,n=n, conf.level=conf.level, alternative=alternative)}
  if(method=="Score"){int=binWilson(Y=Y,n=n, conf.level=conf.level, alternative=alternative)}
  if(method=="AC"){int=binAC(Y=Y,n=n, conf.level=conf.level, alternative=alternative)}
  if(method=="SOC"){int=binSOC(Y=Y,n=n, conf.level=conf.level, alternative=alternative)}
  if(method=="CP"){int=binCP(Y=Y,n=n, conf.level=conf.level, alternative=alternative)}
  if(method=="Blaker"){int=binBlaker(Y=Y,n=n, conf.level=conf.level)}

  if(alternative=="less")
    {CIlength <- int[[2]]-p}

  if(alternative=="greater")
    {CIlength <- p-int[[1]]}

  if(alternative=="two.sided")
    {CIlength <- int[[2]]-int[[1]]}
  CIlength
 }

# Probability of a single event, the binomial density:

 bin.prob <- function(Y,n,p)
  {
   exp( lchoose(n,Y) + Y*log(p) + (n-Y)*log(1-p) )
  }

#  calculate this for all possible events: 

yvec<-0:n

   Lvec<-numeric(length=length(yvec))   
   probvec<-numeric(length=length(yvec))

   for(i in 1:length(yvec))
    {Lvec[i]<-L.Ind.bin(Y=yvec[i], n=n, p=p, conf.level=conf.level, alternative=alternative, method=method)
     probvec[i]<-bin.prob(Y=yvec[i], n=n, p=p)
    }
  expCILength=sum(Lvec * probvec)

# E(X)= sum(Xi * prob(Xi))

out<-list(expCIWidth=expCILength, alternative=alternative, p=p,n=n)

class(out)<-"binWidth"
return(out)
}

