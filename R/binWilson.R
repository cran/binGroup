"binWilson" <-
function(n,Y,conf.level=0.95,alternative="two.sided") {
alpha=1-conf.level
t=Y/n
if(alternative =="two.sided"){
	est.int=(Y+(qnorm(1-alpha/2)^2)/2)/(n+(qnorm(1-alpha/2))^2)
	w.se=((qnorm(1-alpha/2))*sqrt(n*t*(1-t)+(qnorm(1-alpha/2)^2)/4))/(n+qnorm(1-alpha/2)^2)
	KI=c( est.int-w.se, est.int+w.se )
KI}
else{if(alternative=="less"){
	est.int=(Y+(qnorm(1-alpha)^2)/2)/(n+(qnorm(1-alpha))^2)
	w.se=((qnorm(1-alpha))*sqrt(n*t*(1-t)+(qnorm(1-alpha)^2)/4))/(n+qnorm(1-alpha)^2)
	KI=c( 0, est.int+w.se )
} 
else{if(alternative=="greater"){
	est.int=(Y+(qnorm(1-alpha)^2)/2)/(n+(qnorm(1-alpha))^2)
	w.se=((qnorm(1-alpha))*sqrt(n*t*(1-t)+(qnorm(1-alpha)^2)/4))/(n+qnorm(1-alpha)^2)
	KI=c( est.int-w.se , 1 )
}
else{stop("argument alternative misspecified")}}}

conf.int = KI
conf.int
}

