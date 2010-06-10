"bgtCI" <-
function(n, s, y, conf.level=0.95, alternative="two.sided", method="CP")

{
if(length(n)!=1 || (n<1 | abs(round(n)-n) > 1e-07)){stop("number of groups n must be specified as a single integer > 0")}
if(length(s)!=1 || (s<1 | abs(round(s)-s) > 1e-07)){stop("group size s must be specified as a single integer > 0")}
if(length(s)!=1 || (y<0 | abs(round(y)-y) > 1e-07)){stop("observed number of positive groups y must be specified as a single integer>0")}
if(y>n) {stop("number of positive tests y can not be greater than number of groups n")}
if(length(conf.level)!=1 || conf.level<0 || conf.level>1){stop("conf.level must be a positive number between 0 and 1")}

method<-match.arg(method, choices=c("CP","Blaker","AC","Score","Wald","SOC"))
alternative<-match.arg(alternative, choices=c("two.sided","less","greater"))

estimate=1-(1-y/n)^(1/s)

switch(method,
"CP"={conf.int<-bgtCP(n=n, s=s, y=y, conf.level=conf.level, alternative=alternative)},

"Blaker"={
 if(alternative=="less" || alternative=="greater"){
    warning("The Blaker CI is inherently two.sided")
    conf.int<-c(NA, NA)}
 if(alternative=="two.sided")
   {conf.int<-bgtBlaker(n=n, s=s, y=y, conf.level=conf.level)}
 },

"Score"={conf.int<-bgtWilson(n=n, s=s, y=y, conf.level=conf.level, alternative=alternative)},

"AC"={conf.int<-bgtAC(n=n, s=s, y=y, conf.level=conf.level, alternative=alternative)},

"Wald"={conf.int<-bgtWald(n=n, s=s, y=y, conf.level=conf.level, alternative=alternative)},

"SOC"={conf.int<-bgtSOC(n=n, s=s, y=y, conf.level=conf.level, alternative=alternative)}
)

out<-list(conf.int=conf.int,
	estimate=estimate,
	method=method,
	conf.level=conf.level,
	alternative=alternative)
class(out)<-"bgtCI"
out
}

