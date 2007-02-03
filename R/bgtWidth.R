"bgtWidth" <-
function(n, s, p, conf.level=0.95, alternative="two.sided", method="CP")
{
 ns<-n*s
 matnsp<-cbind(n,s,ns,p)

 expCIwidth<-numeric(length=length(matnsp[,1]))

 for (i in 1:length(expCIwidth))
  {
   expCIwidth[i]<-bgtWidthI(n=matnsp[[i,1]], s=matnsp[[i,2]], p=matnsp[[i,4]], conf.level=conf.level, alternative=alternative, method=method)$expCIWidth
  }

 return(as.matrix(cbind(matnsp,expCIwidth)))
}

