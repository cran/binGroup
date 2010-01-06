"bgtAC" <-
function(n, y, s, conf.level=0.95, alternative="two.sided") {
alpha=1-conf.level

est.int=(y+(qnorm(1-alpha/2)^2)/2)/(n+(qnorm(1-alpha/2))^2)
est.int1s=(y+(qnorm(1-alpha)^2)/2)/(n+(qnorm(1-alpha))^2)


if(alternative =="two.sided"){
AC.se=(qnorm(1-alpha/2))*sqrt((est.int*(1-est.int))/(n+(qnorm(1-alpha/2))^2))
KI.int.l=est.int-AC.se
KI.int.u=est.int+AC.se
        if (KI.int.u>1){KI.int.u=1}
        if (KI.int.l<0){KI.int.l=0}
CI=c(1-(1-KI.int.l)^(1/s),1-(1-KI.int.u)^(1/s))
} 

else{if(alternative=="less"){
AC.se=(qnorm(1-alpha))*sqrt((est.int1s*(1-est.int1s))/(n+(qnorm(1-alpha))^2))
KI.int.u=est.int1s+AC.se
        if (KI.int.u>1){KI.int.u=1}
CI=c(0, 1-(1-KI.int.u)^(1/s))
} 

else{if(alternative=="greater"){
AC.se=(qnorm(1-alpha))*sqrt((est.int1s*(1-est.int1s))/(n+(qnorm(1-alpha))^2))
KI.int.l=est.int1s-AC.se
        if (KI.int.l<0){KI.int.l=0}
CI=c(1-(1-KI.int.l)^(1/s),1)
}

else{stop("argument alternative misspecified")}}}

CI
}

