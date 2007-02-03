"binCI" <-
function(n, Y, conf.level=0.95, alternative="two.sided", method="CP")

{
if( n<1 || length(n)!=1){stop("number of groups n must be specified as a single integer > 0")}
if( Y<0 || length(Y)!=1){stop("observed number of positive groups Y must be specified as a single integer>0")}
if(Y>n) {stop("number of positive tests Y can not be greater than number of groups n")}
if(conf.level<0 || conf.level>1 || length(conf.level)!=1){stop("conf.level must be a positive number between 0 and 1, usually 0.95")}
if(method!="CP" && method!="Blaker"&& method!="AC"&& method!="Score"&& method!="Wald"&& method!="SOC"){stop("argument method mis-specified")}
if(alternative!="less" && alternative!="greater"&& alternative!="two.sided"){stop("argument alternative mis-specified")}

estimate=Y/n

 if(method=="CP")
 {
 conf.int=binCP(n=n, Y=Y, conf.level=conf.level, alternative=alternative)
 }

 if(method=="Blaker")
 {
 if(alternative=="less" || alternative=="greater"){cat("the Blaker CI is inherently two.sided","\n")}
 if(alternative=="two.sided")
   {conf.int=binBlaker(n=n, Y=Y, conf.level=conf.level, alternative=alternative)}
 }
 
 if(method=="Score")
 {
 conf.int=binWilson(n=n, Y=Y, conf.level=conf.level, alternative=alternative)
 }

 if(method=="AC")
 {
 conf.int=binAC(n=n, Y=Y, conf.level=conf.level, alternative=alternative)
 }

 if(method=="Wald")
 {
 conf.int=binWald(n=n, Y=Y, conf.level=conf.level, alternative=alternative)
 }

 if(method=="SOC")
 {
 conf.int=binSOC(n=n, Y=Y, conf.level=conf.level, alternative=alternative)
 }

out<-list(conf.int=conf.int,
	estimate=estimate,
	method=method,
	conf.level=conf.level,
	alternative=alternative)
class(out)<-"bgtCI"
out
}

