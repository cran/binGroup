"binAC" <-
function(n, y, conf.level=0.95, alternative="two.sided")
{
alpha=1-conf.level
est=y/n
z1s=qnorm(conf.level)
z2s=qnorm(1-alpha/2)

esti1s=(y+(z1s^2)/2)/(n+z1s^2)
esti2s=(y+(z2s^2)/2)/(n+z2s^2)

ni1s=n+z1s^2
ni2s=n+z2s^2

if(alternative=="two.sided"){

CI=c(esti2s-z2s*sqrt(esti2s*(1-esti2s)/(ni2s)),
     esti2s+z2s*sqrt(esti2s*(1-esti2s)/(ni2s)) )
}
else{if (alternative=="less"){
CI=c( 0 , esti1s+z1s*sqrt(esti1s*(1-esti1s)/(ni1s)) )
}

else{if(alternative=="greater"){
CI=c(esti1s-z1s*sqrt(esti1s*(1-esti1s)/(ni1s)), 1 )
}
else {stop("alternative mis-specified")}}}

CI
}

