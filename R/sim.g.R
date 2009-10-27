sim.g <-
function(gshape=20, gscale=2, beta.par, number.sample, group.size, sens=1, spec=1) {  
  
     covariate<-sort(rgamma(n = number.sample, shape = gshape, scale = gscale))
     prob<-plogis(q = cbind(1, covariate)%*%beta.par)
     ind.binary<-rbinom(n = number.sample, size = 1, prob = prob)
     
     lowest<-floor(number.sample/group.size)
     save1<-matrix(ind.binary[1:(lowest*group.size)], nrow = group.size, ncol = lowest) 
     save.sum<-apply(save1, MARGIN = 2, FUN = sum)
     save.group<-ifelse(save.sum > 0, 1, 0)
     groupresp.t<-ifelse(save.group==1, rbinom(1, 1, sens), 1-rbinom(1, 1, spec))
     groupresp<-rep(groupresp.t, times = 1, each = group.size) 
     group.numb<-rep(1:lowest, times = 1, each = group.size) 

     groupresp.extra<-numeric(0)
     group.numb.extra<-numeric(0)

     if (number.sample/group.size != lowest) {
    
         indice<-(lowest*group.size+1):number.sample
         save.sum.extra<-sum(ind.binary[indice])               
         save.group.extra<-ifelse(save.sum.extra > 0, 1, 0)
         save.group.t.extra<-ifelse(save.group.extra==1, rbinom(1, 1, sens), 1-rbinom(1, 1, spec))
         groupresp.extra<-rep(save.group.t.extra, times = 1, each = length(indice))
         group.numb.extra<-rep(lowest+1, times = 1, each = length(indice)) 
  
     }
     data.frame(groupres = c(groupresp, groupresp.extra), x=covariate, 
                gnum=c(group.numb, group.numb.extra), ind=ind.binary)
}

