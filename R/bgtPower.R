"bgtPower" <-
function(n, s, delta, p.hyp, conf.level=0.95, method="CP", alternative="two.sided")
{

 if( any(n<=3) )
  {stop("the number of groups n allowed in calculations must be integers greater than 1")}
 
 if( any(s<1) ){stop("group size s must be specified as integers > 0")}

 if( length(conf.level)!=1 || conf.level<0 || conf.level>1)
  {stop("conf.level must be a positive number between 0 and 1")}

 if( length(p.hyp)!=1 || p.hyp>1 || p.hyp<0)
  {stop("true proportion p.hyp must be specified as a single number between 0 and 1")}

  method<-match.arg(method, choices=c("CP","Blaker","AC","Score","Wald","SOC"))

  alternative<-match.arg(alternative, choices=c("two.sided","less","greater"))

 if(alternative=="less")
  {if( any( p.hyp-delta < 0) || any(p.hyp-delta > 1) )
   {stop("alternative=less: specify delta as a number between 0 and the threshold p.hyp")}
  }

 if(alternative=="greater")
  {if( any( p.hyp+delta < 0) || any(p.hyp+delta > 1) )
   {stop("alternative=greater: specify delta as a number between the threshold p.hyp and 1")}
  }

 if(alternative=="two.sided")
  {if( any(p.hyp+delta < 0) || any(p.hyp+delta > 1) || any(p.hyp-delta < 0) || any(p.hyp-delta > 1))
   {stop("alternative=two.sided: specify delta as a number between the threshold p.hyp and 1")} 
  }

# calculations:

 matnsp <- cbind(n,s,delta)
 matnsp <- cbind("ns"=matnsp[,1]*matnsp[,2], matnsp)
 power <- numeric(length=nrow(matnsp))
 bias <- numeric(length=nrow(matnsp))


 for( i in 1:length(power))
  {
  temp <- bgtPowerI(n=matnsp[i,2], s=matnsp[i,3], delta=matnsp[i,4], p.hyp=p.hyp, conf.level=conf.level, method=method, alternative=alternative)
  power[i] <- temp$power
  bias[i] <- temp$bias

  }
return(cbind(matnsp, power, bias))

}

