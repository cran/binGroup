

#############################################################################
# This file contains functions to calculate the expected testing expenditure 
# of array testing, and the individual specific measures of accuracy Pooling 
# Sensitivity, Pooling Specificity, Pooling Positive Predicitive Value, and 
# Pooling Negative Predicitive Value.
# Last modified date: 8-15-2011
# Author: Chris McMahan
#############################################################################
#
#Input:
# p: A matrix of probabilities corresponding to the individual's risk
# se: Sensitivity
# sp: Specificity
#
#Output:
# T: Expected testing expenditure of the array
# PSe: A matrix containg each individual's Pooling Sensitivity, corresponds to the input matrix p.
# PSp: A matrix containg each individual's Pooling Specificity, corresponds to the input matrix p.
# PPV: A matrix containg each individual's Pooling Positive Predictive Value, corresponds to the input matrix p.
# NPV: A matrix containg each individual's Pooling Negative Predictive Value, corresponds to the input matrix p.
# (Caution: The function returns the values of PPV and NPV for all individuals in the array even though these measures are
#  diagnostic specific; i.e., PPV (NPV) should only be considered for those individuals who have tested positive (negative),
#  this should be taken into consideration.) 
# 
#
Array.Measures<-function(p,se,sp){
d<-dim(p)
J<-d[1]
K<-d[2]

#######################################################
# Finds probability that each individual tests positive
e.ind<-matrix(-100,nrow=J,ncol=K)
r<-matrix(0:J,ncol=1,nrow=(J+1))
RA<-apply(r,1,rc.pos,p=p,id=1)
c<-matrix(0:K,ncol=1,nrow=(K+1))
CB<-apply(c,1,rc.pos,p=p,id=2)

  for(j in 1:J){
    for(k in 1:K){

      p.temp1<-p
      p.temp1[,k]<-rep(0,J)
      RAJ<-apply(r,1,rc.pos,p=p.temp1,id=1)
      p.c<-prod((1-p[,k]))
      a1<-sum( (1-sp)*(1-se)^r*sp^(J-r)*p.c*RAJ + se*(1-se)^r*sp^(J-r)*((RA-p.c*RAJ)) )

      p.temp1<-p
      p.temp1[j,]<-rep(0,K)
      CBK<-apply(c,1,rc.pos,p=p.temp1,id=2)
      p.r<-prod((1-p[j,]))
      a2<-sum( (1-sp)*(1-se)^c*sp^(K-c)*p.r*CBK + se*(1-se)^c*sp^(K-c)*((CB-p.r*CBK)) )

      a3<- se^2+(1-se-sp)^2*prod(1-p[,k])*prod(1-p[j,])/(1-p[j,k])+(se*(1-sp)-se^2)*(prod(1-p[j,])+prod(1-p[,k]))

      e.ind[j,k]<-(a1+a2+a3)
    }
  }
T<-sum(e.ind)+J+K

###############################################
# Finds Pooling Sensitivity for each individual
pse.ind<-matrix(-100,nrow=J,ncol=K)
c0<-1-(se+(1-se-sp)*apply((1-p),2,prod))
r0<-1-(se+(1-se-sp)*apply((1-p),1,prod))

  for(j in 1:J){
    for(k in 1:K){
      pse.ind[j,k]<-se^3+se^2*(1-se)*(prod(c0[-k])+prod(r0[-j]))
    }
  }

###############################################
# Finds Pooling Specificity for each individual
psp.ind<-matrix(-100,nrow=J,ncol=K)
r<-matrix(0:J,ncol=1,nrow=(J+1))
c<-matrix(0:K,ncol=1,nrow=(K+1))

  for(j in 1:J){
    for(k in 1:K){
      p.temp1<-p
      p.temp1[,k]<-rep(0,J)
      p.temp2<-p
      p.temp2[j,k]<-0

      RAJ<-apply(r,1,rc.pos,p=p.temp1,id=1)
      RAj<-apply(r,1,rc.pos,p=p.temp2,id=1)
      p.c<-prod((1-p.temp2[,k]))
      a1<-sum((1-sp)*(1-se)^r*sp^(J-r)*p.c*RAJ + se*(1-se)^r*sp^(J-r)*((RAj-p.c*RAJ)))

      p.temp1<-p
      p.temp1[j,]<-rep(0,K)
      p.temp2<-p
      p.temp2[j,k]<-0

      CBK<-apply(r,1,rc.pos,p=p.temp1,id=2)
      CBk<-apply(r,1,rc.pos,p=p.temp2,id=2)
      p.r<-prod((1-p.temp2[j,]))
      a2<-sum((1-sp)*(1-se)^c*sp^(K-c)*p.r*CBK + se*(1-se)^c*sp^(K-c)*((CBk-p.r*CBK)))

      a3<-(se+(1-se-sp)*prod((1-p[,k]))/(1-p[j,k]))*(se+(1-se-sp)*prod((1-p[j,]))/(1-p[j,k])) 
      psp.ind[j,k]<-1-(1-sp)*(a1+a2+a3)
    }
  }

#########################################################################
# Finds Pooling Positive (Negative) Predicitive Value for each individual
ppv.ind<-p*pse.ind/(p*pse.ind+(1-p)*(1-psp.ind))
npv.ind<-(1-p)*psp.ind/((1-p)*psp.ind + p*(1-pse.ind))

return(list("T"=T,"PSe"=pse.ind,"PSp"=psp.ind, "PPV"=ppv.ind, "NPV"=npv.ind))
}

##################################################################
##################################################################
# Support Functions needed for Array.Measures()               ####
##################################################################
##################################################################
rc.pos<-function(n,p,id){
d<-dim(p)
M<-d[id]
p.pos<-1-apply((1-p),id,prod)
p.frac<-p.pos/(1-p.pos)
  if(n==0){
  res<-prod(1-p.pos) 
  }
  if(n>0){
    if(n==1){
    res<-sum(p.frac)*prod(1-p.pos)
    }
    if(n>1){
    temp<-cumsum(p.frac[M:n])
    temp<-temp[(M-n+1):1]
      if(n>2){
        for(i in 1:(n-2)){
        temp<-p.frac[(n-i):(M-i)]*temp
        temp<-cumsum(temp[(M-n+1):1])
        temp<-temp[(M-n+1):1]
        }
      }
    temp<-sum(p.frac[1:(M-n+1)]*temp)
    res<-temp*prod(1-p.pos)
    }
  }
return(res)
}


###############################################################################
###############################################################################
###### Spiral or Gradient Array Construction Function                 #########
###############################################################################
###############################################################################
#
# Input:
# prob.vec: A vector of probabilities (of length nr*nc) 
# nr: Number of rows of the array
# nc: Number of columns of the array
# method: Determines whether a spiral (sd) or gradient (gd) is produced
#
# Output:
# array.probs: A matrix of probabilities arranged according to method
Informative.array.prob<-function(prob.vec,nr,nc, method="sd"){

prob.vec<-sort(prob.vec,decreasing=TRUE)
  if(method=="sd"){
    array.probs<-prob.vec[1:2]
    prob.vec<-prob.vec[-(1:2)]
    array.probs<-cbind(array.probs,sort(prob.vec[1:2],decreasing=FALSE))
    prob.vec<-prob.vec[-(1:2)]
     
    max.iter<-max(nr,nc)

      for(i in 1:max.iter){
        if(nrow(array.probs) < nr){
          array.probs<-rbind(array.probs,prob.vec[1:(ncol(array.probs))])
          prob.vec<-prob.vec[-(1:(ncol(array.probs)))]
        }
        if(ncol(array.probs) < nc){
          array.probs<-cbind(array.probs,sort(prob.vec[1:(nrow(array.probs))],decreasing=FALSE))
          prob.vec<-prob.vec[-(1:(nrow(array.probs)))]
        }
      }  
    }

  if(method=="gd"){
    array.probs<-matrix(prob.vec,ncol=max(nr,nc), nrow=min(nr,nc), byrow=FALSE)
  }
return(array.probs)
}





######################################################################################################
# Example 1: This is an example of finding the expected testing expenditure and measures of accuracy #
#            for a homogeneous population, note the results are the same as those found using the    #
#            equations derived in Kim et al. 2007.                                                   #
######################################################################################################


se=0.90     # Testing Sensitivity
sp=0.99     # Testing Specificity
p=0.045     # Population prevalence
J=10        # The number of rows
K=10        # The number of columns

pmat<-matrix(p,nrow=J,ncol=K) # Builds matrix of probabilities 
Array.Measures(p=pmat,se=se,sp=sp)



######################################################################################################
# Example 2: This is an example of finding the expected testing expenditure and measures of accuracy #
#            for a heterogeneous population.                                                         #
######################################################################################################


se=0.90     # Testing Sensitivity
sp=0.99     # Testing Specificity
p=0.045     # Population mean prevalence
alpha<-1    # Controls the amount of heterogeneity in the population; i.e., heterogeneity increases as alpha decreases and vice versa  
J=10        # The number of rows
K=10        # The number of columns

p<-rbeta(J*K,alpha,((alpha-alpha*p)/p))

p.ra<-matrix(p,nrow=J,ncol=K)                                      # Randomly filling the array
p.ga<-Informative.array.prob(prob.vec=p,nr=J,nc=K, method="gd")    # Filling the array according to the gradient design
p.sa<-Informative.array.prob(prob.vec=p,nr=J,nc=K, method="sd")    # Filling the array according to the spiral design

Array.Measures(p=p.ra,se=se,sp=sp)
Array.Measures(p=p.ga,se=se,sp=sp)
Array.Measures(p=p.sa,se=se,sp=sp)











