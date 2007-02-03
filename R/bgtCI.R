"bgtCI" <-
function(n, s, Y, conf.level=0.95, alternative="two.sided", method="CP")

{
if( n<1 || length(n)!=1){stop("number of groups n must be specified as a single integer > 0")}
if( s<1 || length(s)!=1){stop("group size s must be specified as a single integer > 0")}
if( Y<0 || length(Y)!=1){stop("observed number of positive groups Y must be specified as a single integer>0")}
if(Y>n) {stop("number of positive tests Y can not be greater than number of groups n")}
if(conf.level<0 || conf.level>1 || length(conf.level)!=1){stop("conf.level must be a positive number between 0 and 1, usually 0.95")}
if(method!="CP" && method!="Blaker"&& method!="AC"&& method!="Score"&& method!="Wald"&& method!="SOC"){stop("argument method mis-specified")}
if(alternative!="less" && alternative!="greater"&& alternative!="two.sided"){stop("argument alternative mis-specified")}

estimate=1-(1-Y/n)^(1/s)

 if(method=="CP")
 {
 conf.int=bgtCP(n=n, s=s, Y=Y, conf.level=conf.level, alternative=alternative)
 }

 if(method=="Blaker")
 {
 if(alternative=="less" || alternative=="greater"){cat("the Blaker CI is inherently two.sided","\n")}
 if(alternative=="two.sided")
   {conf.int=bgtBlaker(n=n, s=s, Y=Y, conf.level=conf.level, alternative=alternative)}
 }
 
 if(method=="Score")
 {
 conf.int=bgtWilson(n=n, s=s, Y=Y, conf.level=conf.level, alternative=alternative)
 }

 if(method=="AC")
 {
 conf.int=bgtAC(n=n, s=s, Y=Y, conf.level=conf.level, alternative=alternative)
 }

 if(method=="Wald")
 {
 conf.int=bgtWald(n=n, s=s, Y=Y, conf.level=conf.level, alternative=alternative)
 }

 if(method=="SOC")
 {
 conf.int=bgtSOC(n=n, s=s, Y=Y, conf.level=conf.level, alternative=alternative)
 }

out<-list(conf.int=conf.int,
	estimate=estimate,
	method=method,
	conf.level=conf.level,
	alternative=alternative)
class(out)<-"bgtCI"
out
}

