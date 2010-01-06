"bgtWidth" <-
function(n, s, p, conf.level=0.95, alternative="two.sided", method="CP")
{

 if( any(n<=3) )
  {stop("the number of groups n allowed in calculations must be integers greater than 1")}
 
 if( any(s<1) ){stop("group size s must be specified as integers > 0")}

 if( length(conf.level)!=1 || conf.level<0 || conf.level>1)
  {stop("conf.level must be a positive number between 0 and 1")}

 if( length(p)!=1 || p>1 || p<0)
  {stop("true proportion p must be specified as a single number between 0 and 1")}

  method<-match.arg(method, choices=c("CP","Blaker","AC","Score","Wald","SOC"))

  alternative<-match.arg(alternative, choices=c("two.sided","less","greater"))

# calculations:

 matnsp <- cbind(n,s,p)
 matnsp <- cbind("ns"=matnsp[,1]*matnsp[,2], matnsp)
 power <- numeric(length=nrow(matnsp))
 bias <- numeric(length=nrow(matnsp))

 expCIwidth<-numeric(length=nrow(matnsp))

 for (i in 1:length(expCIwidth))
  {
   expCIwidth[i]<-bgtWidthI(n=matnsp[[i,2]], s=matnsp[[i,3]], p=matnsp[[i,4]], conf.level=conf.level, alternative=alternative, method=method)$expCIWidth
  }

 return(as.matrix(cbind(matnsp,expCIwidth)))
}

