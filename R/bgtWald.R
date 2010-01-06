"bgtWald" <-
function (n, y, s, conf.level=0.95, alternative="two.sided")

{
if(y>n) {stop("number of positive tests y can not be greater than number of tests n")}
th=y/n
esti=1-(1-th)^(1/s)
var.esti=(1-(1-esti)^s)/(n*(s^2)*(1-esti)^(s-2))
alpha=1-conf.level

if(alternative=="two.sided"){
    snquant=qnorm(p=1-alpha/2,mean=0,sd=1,lower.tail=TRUE)
    CI=c(esti-snquant*sqrt(var.esti),esti+snquant*sqrt(var.esti))
}
else{if (alternative=="less"){
    snquant=qnorm(p=1-alpha,mean=0,sd=1,lower.tail=TRUE)
    CI=c(0 ,esti+snquant*sqrt(var.esti))
}
else {if (alternative=="greater"){
    snquant=qnorm(p=1-alpha,mean=0,sd=1,lower.tail=TRUE)
    CI=c(esti-snquant*sqrt(var.esti), 1)
}
else {stop("argument alternative mis-specified")}}}
CI
}

