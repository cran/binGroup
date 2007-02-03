"bgtPowerI" <-
function(n, s, delta, p.hyp, conf.level, method, alternative){


# 1) P.Ind returns TRUE in case that the CI doesnot contain p.hyp for a certain event Y=y


P.Ind <- function(n,Y,s,p.hyp,conf.level,method,alternative){

 if(method=="Score"){
  KI.Wilson <- bgtWilson(n=n,Y=Y,s=s,conf.level=conf.level, alternative=alternative) 

   (KI.Wilson[[1]]>=p.hyp||KI.Wilson[[2]]<=p.hyp)
  }

 else{if(method=="AC"){
  KI.AC<-bgtAC(n=n,Y=Y,s=s,conf.level=conf.level, alternative=alternative) 

  (KI.AC[[1]]>=p.hyp||KI.AC[[2]]<=p.hyp)
  }

 else{if(method=="Wald"){
  KI.Wald<-bgtWald(n=n,Y=Y,s=s,conf.level=conf.level, alternative=alternative)

  (KI.Wald[[1]]>=p.hyp||KI.Wald[[2]]<=p.hyp)
  }

 else{if(method=="CP"){
  KI.CP<-bgtCP(n=n,Y=Y,s=s,conf.level=conf.level, alternative=alternative) 

  (KI.CP[[1]]>=p.hyp||KI.CP[[2]]<=p.hyp)
  }

 else{if(method=="SOC"){
  KI.SOC<-bgtSOC(n=n,Y=Y,s=s,conf.level=conf.level, alternative=alternative)

  (KI.SOC[[1]]>=p.hyp||KI.SOC[[2]]<=p.hyp)
 }

 else{if(method=="Blaker"){
  KI.Bl<-bgtBlaker(n=n,Y=Y,s=s,conf.level=conf.level) 

   if(alternative=="two.sided")
   {dec<-(KI.Bl[[1]]>=p.hyp||KI.Bl[[2]]<=p.hyp)}

   if(alternative=="less")
   {dec<-(KI.Bl[[2]]<=p.hyp)}

   if(alternative=="greater")
   {dec<-(KI.Bl[[1]]>=p.hyp)}
 dec
 }

 else{stop("argument method mis-specified")}}}}}}
}

# end of P.Ind

# # #

# 2) Probability of a certain event Y=y:

 bgt.prob<-function(n,Y,s,p.tr)
  {
  theta<-1-(1-p.tr)^s
  dbinom(x=Y,size=n, prob=theta)
  }

# # # 
#  3) sum( P.Ind(Y=y)*Prob(Y=y) ) for all realizations of y:

if(alternative=="less" || alternative=="greater")
 {
 
 if(alternative=="less"){p.tr = p.hyp - delta}
 if(alternative=="greater"){p.tr = p.hyp + delta}

 yvec<-0:n

 probvec<-numeric(length=length(yvec))
 powvec<-numeric(length=length(yvec))
 expvec<-numeric(length=length(yvec))

 for(i in 1:length(yvec))
  {
  probvec[i] <- bgt.prob(n=n,Y=yvec[i],s=s,p.tr=p.tr)
  powvec[i] <- P.Ind(n=n,Y=yvec[i],s=s,p.hyp=p.hyp,conf.level=conf.level,method=method, alternative = alternative)
  expvec[i] <- (1-(1-yvec[i]/n)^(1/s))
  }
 powex<-sum(powvec * probvec) 
 expex<-sum(expvec * probvec)
 bias<-expex-p.tr

 out<-list(power=powex, bias=bias, p.tr=p.tr)

 }

if(alternative=="two.sided")
 {
 p.trl = p.hyp - delta
 p.trg = p.hyp + delta

 yvec<-0:n
 powvec<-numeric(length=length(yvec))
 expvec<-numeric(length=length(yvec))
 probvecl<-numeric(length=length(yvec))
 probvecg<-numeric(length=length(yvec))



 for(i in 1:length(yvec))
 {
  probvecl[i] <- bgt.prob(n=n,Y=yvec[i],s=s,p.tr=p.trl)
  probvecg[i] <- bgt.prob(n=n,Y=yvec[i],s=s,p.tr=p.trg)

  powvec[i] <- P.Ind(n=n,Y=yvec[i],s=s,p.hyp=p.hyp,conf.level=conf.level,method=method, alternative = alternative)

  expvec[i] <- (1-(1-yvec[i]/n)^(1/s))
 }

 powexl<-sum(powvec * probvecl) 
 expexl<-sum(expvec * probvecl)
 biasl<-expexl-p.trl

 powexg<-sum(powvec * probvecg) 
 expexg<-sum(expvec * probvecg)
 biasg<-expexg-p.trg

 out<-list(power=min(powexl,powexg),bias=max(biasl,biasg), p.tr=c(p.trl,p.trg))

 }

class(out)<-"bgtpower"
out

}

