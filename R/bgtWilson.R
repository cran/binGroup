"bgtWilson" <-
function(n, y, s, conf.level=0.95, alternative="two.sided")
{ 
alpha=1-conf.level 
th=y/n
est.int=(y+(qnorm(1-alpha/2)^2)/2)/(n+(qnorm(1-alpha/2))^2)
est.int1s=(y+(qnorm(1-alpha)^2)/2)/(n+(qnorm(1-alpha))^2)

if(alternative =="two.sided"){
    w.se=((qnorm(1-alpha/2))*sqrt(n*th*(1-th)+(qnorm(1-alpha/2)^2)/4))/(n+qnorm(1-alpha/2)^2)
    KI.int.l=est.int-w.se
    KI.int.u=est.int+w.se
        if (KI.int.u>1){KI.int.u=1}
        if (KI.int.l<0){KI.int.l=0}
    KI=c( 1-(1-KI.int.l)^(1/s), 1-(1-KI.int.u)^(1/s) )
}

else{if(alternative=="less"){
    w.se=((qnorm(1-alpha))*sqrt(n*th*(1-th)+(qnorm(1-alpha)^2)/4))/(n+qnorm(1-alpha)^2)
    KI.int.u=est.int1s+w.se
        if (KI.int.u>1){KI.int.u=1}
    KI=c( 0, 1-(1-KI.int.u)^(1/s) )
}

else{if(alternative=="greater"){
    w.se=((qnorm(1-alpha))*sqrt(n*th*(1-th)+(qnorm(1-alpha)^2)/4))/(n+qnorm(1-alpha)^2)
    KI.int.l=est.int1s-w.se
        if (KI.int.l<0){KI.int.l=0}
    KI=c( 1-(1-KI.int.l)^(1/s), 1 )
}

else{stop("argument alternative misspecified")}}}

KI

}

