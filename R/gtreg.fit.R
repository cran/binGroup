gtreg.fit <-
function(resp, cova, groupn, sens, spec, linkf, 
                  optim.meth="Nelder-Mead") {

    z<-tapply(resp, groupn, tail, n=1)
    K<-ncol(cova)
    
    if (K==1) cova.mean<-tapply(cova, groupn, mean) else {
         temp<-by(cova, groupn, mean)
         cova.mean<-do.call(rbind,temp)
    }
       
    mod.fit<-lm.fit(x=cova.mean, y=z) #use group mean age to predict Z
    beta.group<-mod.fit$coefficients

    logL<-function(beta, x) {

          p.ijk<-switch(linkf, logit=plogis(x%*%beta), probit=pnorm(x%*%beta),
                cloglog=1-exp(-exp(x%*%beta)))

          prod.p.ijk<-tapply(1-p.ijk, groupn, prod)

          -sum(z*log(sens+(1-sens-spec)*prod.p.ijk) + (1-z)*log(1-sens-(1-sens-spec)*prod.p.ijk))

    }
    
    group.mod.fit<-optim(par = beta.group, fn = logL, method = optim.meth,  
               x=cova, control = list(trace = 0), hessian = TRUE)

    mod.fit.Ho<-lm(formula = z ~ 1)
    mod.fit0<-optim(par = mod.fit.Ho$coefficients, fn = logL, method = "BFGS",  
               x=matrix(1, nrow=length(resp)), control = list(trace = 0, maxit = 1000))
    nulld<-2*mod.fit0$value
    residd<-2*group.mod.fit$value

    lin.pred<-cova%*%group.mod.fit$par

    p.ijk<-switch(linkf, logit=plogis(lin.pred), probit=pnorm(lin.pred),
         cloglog=1-exp(-exp(lin.pred)))
    
    prod.p.ijk<-tapply(1-p.ijk, groupn, prod)
    
    zhat<-sens+(1-sens-spec)*prod.p.ijk
    residual<-z-zhat

    aic<-residd+2*K

    if (group.mod.fit$convergence==0) counts<-group.mod.fit$counts[[1]]
        else counts<-"did not converge"
    
    list(coefficients=group.mod.fit$par, hessian=group.mod.fit$hessian, fitted.group.values=zhat,
           deviance=residd, df.residual=length(resp)-K, null.deviance=nulld, df.null=length(resp)-1,
           aic=aic, counts=counts, residuals=residual, z=z)

}

