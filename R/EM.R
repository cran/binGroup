EM <-
function(resp, cova, groupn, sens, spec, linkf, EM.maxiter=1000)  {

    z<-tapply(resp, groupn, tail, n=1)
    K<-ncol(cova)
    
    if (K==1) cova.mean<-tapply(cova, groupn, mean) else {
         temp<-by(cova, groupn, mean)
         cova.mean<-do.call(rbind,temp)
    }
    
    sam<-length(resp)
    vec<-1:sam
    
    group.sizes<-tapply(resp, groupn, length)

    diff<-1
    counts<-1
    
    mod.fit<-lm.fit(x=cova.mean, y=z) #use group mean age to predict Z
    beta.old<-mod.fit$coefficients

    while(diff > 0.001 & counts <= EM.maxiter) { 
  
        p.ijk<-switch(linkf, logit=plogis(cova%*%beta.old), probit=pnorm(cova%*%beta.old),
             cloglog=1-exp(-exp(cova%*%beta.old)))
        
        prod.p.ijk<-tapply(1-p.ijk, groupn, prod)
        
        #E-Step: Find E(Y_ik | Z_k)                                                                     
        den<-rep((1-spec)*prod.p.ijk+sens*(1-prod.p.ijk), group.sizes)
        den2<-rep(spec*prod.p.ijk+(1-sens)*(1-prod.p.ijk), group.sizes)
        expect<-rep(NA, times=sam)
      
        for (counter in vec) {
            if (resp[counter] == 0) expect[counter] <- (1-sens)*p.ijk[counter]/den2[counter]
                else expect[counter]<- sens*p.ijk[counter]/den[counter]
        }
    
    
        #M-step: Fit model  
       
        binlogL<-function(beta) {
         
            pij<-switch(linkf, logit=plogis(cova%*%beta), probit=pnorm(cova%*%beta),
                  cloglog=1-exp(-exp(cova%*%beta)))
            
            -sum(expect*log(pij) + (1-expect)*log(1-pij))
        }
       
        mod.fit<-optim(par = beta.old, fn = binlogL)
        diff<-max(abs((beta.old - mod.fit$par)/beta.old))
        beta.old<-mod.fit$par
        counts<-counts+1
   
    }
   
   
    #calculate Hessian matrix
   
    xib<-cova%*%beta.old
    erf<-2*p.ijk-1
    pt1<-switch(linkf, logit=-exp(xib)/(1+exp(xib))^2, 
          probit=sqrt(2)*xib*exp(-xib^2/2)/(sqrt(pi)*(1-erf))-2*exp(-xib^2)/(pi*(1-erf)^2),
          cloglog=-exp(xib))
    pt2<-switch(linkf, logit=0,
          probit=(8*exp(-xib^2/2)*erf+2*xib*sqrt(2*pi)*erf^2-2*xib*sqrt(2*pi))
                *exp(-xib^2/2)/((1+erf)^2*pi*(1-erf)^2),
          cloglog=-(exp(xib-exp(xib))+exp(2*xib-exp(xib))-exp(xib))/(exp(-exp(xib))-1)^2)
    nm<-pt1+expect*pt2
    sign1<-as.vector(sign(nm))
    nn<-as.vector(sqrt(abs(nm)))
    x2<-cova*nn
    m<-(t(x2)%*%(sign1*x2))
   
    b<-array(NA,c(K,K,sum(group.sizes^2)))
    p<-1

    for(i in vec)
       for(j in vec[groupn==groupn[i]]) {
           wii<-ifelse(i==j, expect[i]-expect[i]^2, expect[i]*(p.ijk[j]-expect[j]))
           coe<-switch(linkf, logit=1, 
               probit=8*exp(-(xib[i]^2+xib[j]^2)/2)/((1-erf[i]^2)*(1-erf[j]^2)*pi),
               cloglog=exp(xib[i]+xib[j])/((exp(-exp(xib[i]))-1)*(exp(-exp(xib[j]))-1)))
           b[,,p]<-wii*coe*cova[i,]%*%t(cova[j,])
           p<-p+1 
       }
   
    m1<-apply(b, c(1,2), sum)
    H<--(m+m1)


    #residuals
    zhat<-sens+(1-sens-spec)*prod.p.ijk
    residual<-z-zhat

    #residual deviance
    residd<--2*sum(z*log(zhat) + (1-z)*log(1-zhat))
   
    logL0<-function(beta) {

          inter<-rep(beta, sam)
          p.ijk<-switch(linkf, logit=plogis(inter), probit=pnorm(inter),
                cloglog=1-exp(-exp(inter)))

          prod.p.ijk<-tapply(1-p.ijk, groupn, prod)

          -sum(z*log(sens+(1-sens-spec)*prod.p.ijk) + (1-z)*log(1-sens-(1-sens-spec)*prod.p.ijk))

    }
    
    #null deviances
    mod.fit.Ho<-lm(formula = z ~ 1)
    mod.fit0<-optim(par = mod.fit.Ho$coefficients, fn = logL0, method = "BFGS",  
                  control = list(trace = 0, maxit = 1000))
    nulld<-2*mod.fit0$value
   
    aic<-residd+2*K
   
    if (counts == EM.maxiter) counts<-"Xie did not converge"
   
    list(coefficients=beta.old, hessian=H, fitted.group.values=zhat, deviance=residd, 
          df.residual=sam-K, null.deviance=nulld, df.null=sam-1,
          aic=aic, counts=counts, residuals=residual, z=z)

}

