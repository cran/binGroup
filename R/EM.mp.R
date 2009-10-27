EM.mp <-
function(col.resp, row.resp, cova, coln, rown, sqn, ret, sens, spec, 
        linkf, n.gibbs=1000, burn.num=20, tol=0.005, EM.maxiter=500, init.par=NULL, sens.ind=NULL, spec.ind=NULL)  {
     
     start.time<-proc.time()
     if (is.null(sens.ind)) sens.ind<-sens
     if (is.null(spec.ind)) spec.ind<-spec
     
     len<-max(sqn)
     diff<-1
     count<-1
     sam<-length(sqn)
     
     col.groupn<-coln[sqn==1]
     if (len>1) {
        for (i in 2:len)  { 
            temp<-max(col.groupn)+coln[sqn==i]
            col.groupn<-c(col.groupn, temp)
        }
     }
     
     if (is.null(init.par)) {
        mod.fit<-try(gtreg.fit(col.resp, cova, col.groupn, sens, spec, linkf))
        if (class(mod.fit)=="try-error") {
            row.groupn<-rown[sqn==1]
            if (len>1) {
                for (i in 2:len)  { 
                     temp<-max(row.groupn)+rown[sqn==i]
                     row.groupn<-c(row.groupn, temp)
                }
            }
            mod.fit<-gtreg.fit(row.resp, cova, row.groupn, sens, spec, linkf)        
        }    
        beta.old<-mod.fit$coefficients
     }       
     else {
        beta.old<-init.par
        names(beta.old)<-dimnames(cova)[[2]]
     }
      
     
     while(diff > tol && count <= EM.maxiter) { 
  
           p.ijk.all<-switch(linkf, logit=plogis(cova%*%beta.old), probit=pnorm(cova%*%beta.old),
                      cloglog=1-exp(-exp(cova%*%beta.old)))
           expect.all<-numeric(0)
           
           mat2<-index<-0
           xib<-cova%*%beta.old
           erf<-2*p.ijk.all-1
           
           for (arrayn in 1:len) { #arrayn is the array number
            
               index.r<-index.c<-vector("logical", length=sam)
               for (i in 1:sam) {
                    if (rown[i]==1 && sqn[i]==arrayn) index.c[i]<-TRUE 
                        else index.c[i]<-FALSE
                    if (coln[i]==1 && sqn[i]==arrayn) index.r[i]<-TRUE
                        else index.r[i]<-FALSE
               }
               n.row<-max(rown[index.r])  #number of rows in current array
               n.col<-max(coln[index.c])
               rowresp<-row.resp[index.r] #row.resp for current array
               colresp<-col.resp[index.c]
               
               index<-max(index)+1:(n.row*n.col)
               
               if (!is.null(ret)) {
                  re.ind<-na.omit(cbind(coln[sqn==arrayn], rown[sqn==arrayn], 
                              ret[sqn==arrayn]))
                  re<-ifelse(re.ind[,3]==1, sens.ind, 1-sens.ind)   #  |assume yij=1
                  re1<-ifelse(re.ind[,3]==0, spec.ind, 1-spec.ind)  #  |assume yij=0
               } 
               
               p.ijk<-matrix(p.ijk.all[sqn==arrayn], nrow=n.row)
               a<-ifelse(rowresp==1, sens, 1-sens) #vector operation, sensitivity or not
               b<-ifelse(colresp==1, sens, 1-sens)
               a1<-ifelse(rowresp==0, spec, 1-spec)
               b1<-ifelse(colresp==0, spec, 1-spec)
      
               mat<-array(NA,c(n.row, n.col, n.gibbs))
               y<-matrix(0, nrow=n.row, ncol=n.col) #individual responses
               
               for (k in 1:(n.gibbs+burn.num)) { #the burn-in period
                            l<-1
                            for (j in 1:n.col)   #need to repeat this n.gibbs times
                                for (i in 1:n.row) {
                                    num<-a[i]*b[j]*p.ijk[i,j]
                                    den.r<-ifelse(sum(y[i,])-y[i,j]>0, a[i], a1[i])  #p(r_i|real r_i)
                                    den.c<-ifelse(sum(y[,j])-y[i,j]>0, b[j], b1[j])  #p(c_i|real c_i)
                                    den2<-den.r*den.c*(1-p.ijk[i,j])
                                    if (!is.null(ret)) {
                                       if (l<=length(re) && j==re.ind[l,1] && i==re.ind[l,2]) {
                                           num<-num*re[l]
                                           den2<-den2*re1[l]
                                           l<-l+1
                                       }
                                    } 
                                    den<-num+den2              #conditional prob
                                    if (den!=0) {
                                       cond.p<-num/den
                                       y[i,j]<-rbinom(1,1,cond.p)
                                    }
                                    else y[i,j]<-0
                                }
     
                            if (k>burn.num) {
                            
                                mat[,,k-burn.num]<-y
                                vec<-as.vector(y)
                                
                                for(i1 in index[vec==1])
                                      for(j1 in index[vec==1]) {
                                          bq<-switch(linkf, logit=1, 
                                                 probit=8*exp(-(xib[i1]^2+xib[j1]^2)/2)
                                                 /((1-erf[i1]^2)*(1-erf[j1]^2)*pi),
                                                 cloglog=exp(xib[i1]+xib[j1])/
                                                 ((exp(-exp(xib[i1]))-1)*(exp(-exp(xib[j1]))-1)))*
                                                 cova[i1,]%*%t(cova[j1,])
                                          mat2<-mat2+bq 
                                      }
                            }

               }            
               expect.m<-apply(mat, c(1,2), mean)
               expect<-as.vector(expect.m)
               expect.all<-c(expect.all, expect)
               
           }  #end of array

           #M-step: Fit model  
           binlogL<-function(beta) {
    
               pij<-switch(linkf, logit=plogis(cova%*%beta), probit=pnorm(cova%*%beta),
                     cloglog=1-exp(-exp(cova%*%beta)))
    
               -sum(expect.all*log(pij) + (1-expect.all)*log(1-pij))
           }
           
           mod.fit<-optim(par = beta.old, fn = binlogL)
           diff<-max(abs((beta.old - mod.fit$par)/beta.old)) #might be too harsh
           beta.old<-mod.fit$par
           count<-count+1
                
           cat("beta is:",beta.old,"\tdiff is:",diff,"\n")

     }  #end of beta
     
     #calculate Hessian matrix
     
     index<-0        
     first<-mat2/n.gibbs #p*p matrix
     second<-0
     for (arrayn in 1:len) {     
          
          n.row<-max(rown[sqn==arrayn])
          n.col<-max(coln[sqn==arrayn])        
          index<-max(index)+1:(n.row*n.col)
          expect<-expect.all[index]
          for(i1 in index[expect!=0])
              for(j1 in index[expect!=0]) {
                  coe<-switch(linkf, logit=1, 
                              probit=8*exp(-(xib[i1]^2+xib[j1]^2)/2)
                              /((1-erf[i1]^2)*(1-erf[j1]^2)*pi),
                              cloglog=exp(xib[i1]+xib[j1])/
                              ((exp(-exp(xib[i1]))-1)*(exp(-exp(xib[j1]))-1)))
                  tim<-expect.all[i1]*expect.all[j1]*coe*cova[i1,]%*%t(cova[j1,])
                  second<-second+tim 
              }               
     }         
     p2.array<-first-second

     pt1<-switch(linkf, logit=-exp(xib)/(1+exp(xib))^2, 
            probit=sqrt(2)*xib*exp(-xib^2/2)/(sqrt(pi)*(1-erf))-2*exp(-xib^2)/(pi*(1-erf)^2),
            cloglog=-exp(xib))
     pt2<-switch(linkf, logit=0,
            probit=(8*exp(-xib^2/2)*erf+2*xib*sqrt(2*pi)*erf^2-2*xib*sqrt(2*pi))
                *exp(-xib^2/2)/((1+erf)^2*pi*(1-erf)^2),
            cloglog=-(exp(xib-exp(xib))+exp(2*xib-exp(xib))-exp(xib))/(exp(-exp(xib))-1)^2)
     nm<-pt1+expect.all*pt2
     sign1<-as.vector(sign(nm))
     nn<-as.vector(sqrt(abs(nm)))
     x2<-cova*nn
     m<-(t(x2)%*%(sign1*x2))
     H<--(m+p2.array)

     end.time<-proc.time()

     save.time<-end.time-start.time

     cat("\n Number of minutes running:", save.time[3]/60, "\n \n")
     
     list(coefficients=beta.old, hessian=H, df.residual=sam-ncol(cova), df.null=sam-1, 
              Gibbs.sample.size=n.gibbs, counts=count)

}

