"binWald" <-
function(n, y, conf.level=0.95, alternative="two.sided")
{
alpha=1-conf.level
est=y/n
z1s=qnorm(conf.level)
z2s=qnorm(1-alpha/2)

if(alternative=="two.sided"){
KI=c(est-z2s*sqrt(est*(1-est)/(n)),
     est+z2s*sqrt(est*(1-est)/(n)) )
}
else{if (alternative=="less"){
KI=c( 0 , est+z1s*sqrt(est*(1-est)/(n)) )
}
else{if(alternative=="greater"){
KI=c(est-z1s*sqrt(est*(1-est)/(n)), 1 )
}
else {stop("alternative mis-specified")}}}
KI
}

