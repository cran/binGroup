"bgtPower" <-
function(n, s, delta, p.hyp, conf.level=0.95, method="CP", alternative="two.sided")
{

 if( any(n<=3) )
  {stop("the number of groups n allowed in calculations must be integers greater than 1")}
 
 if( any(s<1) ){stop("group size s must be specified as integers > 0")}

 if(conf.level<0 || conf.level>1 || length(conf.level)!=1)
  {stop("conf.level must be a positive number between 0 and 1, usually 0.95")}

 if(method!="CP" && method!="Blaker"&& method!="AC"&& method!="Score"&& method!="Wald"&& method!="SOC")
  {stop("argument method mis-specified")}

 if(alternative!="less" && alternative!="greater"&& alternative!="two.sided")
  {stop("argument alternative mis-specified")}

 if( p.hyp>1 || p.hyp<0 || length(p.hyp)!=1)
  {stop("true proportion p.hyp must be specified as a single number between 0 and 1")}

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

 ns <- n*s
 matnsp <- cbind(n,s,ns,delta)
 power <- numeric(length=length(matnsp[,1]))
 bias <- numeric(length=length(matnsp[,1]))


 for( i in 1:length(power))
  {
  temp <- bgtPowerI(n=matnsp[i,1], s=matnsp[i,2], delta=matnsp[i,4], p.hyp=p.hyp, conf.level=conf.level, method=method, alternative=alternative)
  power[i] <- temp$power
  bias[i] <- temp$bias

  }
return(cbind(matnsp, power, bias))

}

