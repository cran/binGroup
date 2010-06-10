"binCI" <-
function(n, y, conf.level=0.95, alternative="two.sided", method="CP")

{
if(length(n)!=1 || (n<1 | abs(round(n)-n) > 1e-07)){stop("number of groups n must be specified as a single integer > 0")}
if(length(y)!=1 || (y<0 | abs(round(y)-y) > 1e-07)){stop("observed number of positive groups y must be specified as a single integer>0")}
if(y>n) {stop("number of positive tests y can not be greater than number of groups n")}
if( length(conf.level)!=1 || conf.level<0 || conf.level>1){stop("conf.level must be a positive number between 0 and 1, usually 0.95")}

method<-match.arg(method, choices=c("CP","Blaker","AC","Score","Wald","SOC"))
alternative<-match.arg(alternative, choices=c("two.sided","less","greater"))

estimate=y/n

switch(method,
"CP"={conf.int<-binCP(n=n, y=y, conf.level=conf.level, alternative=alternative)},

"Blaker"={
 if(alternative=="less" || alternative=="greater"){
  warning("The Blaker CI is inherently two-sided.")
  conf.int<-c(NA, NA)
  }
 if(alternative=="two.sided")
  {conf.int<-binBlaker(n=n, y=y, conf.level=conf.level)}
 },
 
"Score"={conf.int<-binWilson(n=n, y=y, conf.level=conf.level, alternative=alternative)},

"AC"={conf.int<-binAC(n=n, y=y, conf.level=conf.level, alternative=alternative)},

"Wald"={conf.int<-binWald(n=n, y=y, conf.level=conf.level, alternative=alternative)},

"SOC"={conf.int<-binSOC(n=n, y=y, conf.level=conf.level, alternative=alternative)}
)

out<-list(conf.int=conf.int,
	estimate=estimate,
	method=method,
	conf.level=conf.level,
	alternative=alternative)
class(out)<-"binCI"
out
}

