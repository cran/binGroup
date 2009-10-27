sim.mp <-
function(gshape=20, gscale=2, beta.par, rown, coln, sens=1, spec=1, sens.ind=NULL, spec.ind=NULL) {

     if (is.null(sens.ind)) sens.ind<-sens
     if (is.null(spec.ind)) spec.ind<-spec

     number.sample<-sum(coln*rown)
     len<-length(rown)
     cova<-rgamma(n = number.sample, shape = gshape, scale = gscale)
     pij<-as.vector(plogis(q = cbind(1, cova)%*%beta.par))
     ind.binary<-rbinom(n = number.sample, size = 1, prob = pij)
     individual<-col.groupn<-row.groupn<-numeric(0)
     rep.row.resp<-rep.col.resp<-numeric(0)
     ret<-rep(NA, number.sample)

     for (i in 1:len) {
         if (i>1) index<-seq(max(index)+1, length=(rown*coln)[i]) else index<-1:(rown*coln)[1]
         ind<-matrix(ind.binary[index], nrow=rown[i])
         col.resp<-apply(ind, MARGIN = 2, FUN = sum)
         col.resp<-ifelse(col.resp>0,1,0)
         col.terror<-ifelse(col.resp==1, rbinom(1, 1, sens), rbinom(1, 1, 1-spec))

         row.resp<-apply(ind, MARGIN = 1, FUN = sum)
         row.resp<-ifelse(row.resp>0,1,0)
         row.terror<-ifelse(row.resp==1, rbinom(1, 1, sens), rbinom(1, 1, 1-spec))

         temp.c<-rep(1:coln[i], each=rown[i])  #col group number
         col.groupn<-c(col.groupn, temp.c)

         temp.r<-rep(1:rown[i], coln[i])  #col group number
         row.groupn<-c(row.groupn, temp.r)

         temp2.c<-rep(col.terror, each=rown[i])  #repeated col group responses
         rep.col.resp<-c(rep.col.resp, temp2.c)

         temp2.r<-rep(row.terror, coln[i])  #repeated col group responses
         rep.row.resp<-c(rep.row.resp, temp2.r)

         if (all(rep.row.resp==0)) {
             for (i in index) {
                  if (rep.col.resp[i]==1) ret[i]<-ifelse(ind.binary[i]==1,
                                                      rbinom(1, 1, sens.ind),
                                                      rbinom(1, 1, 1-spec.ind))
             }
         }

         if (all(rep.col.resp==0)) {
             for (i in index) {
                  if (rep.row.resp[i]==1) ret[i]<-ifelse(ind.binary[i]==1,
                                                      rbinom(1, 1, sens.ind),
                                                      rbinom(1, 1, 1-spec.ind))
             }
         }

         individual<-c(individual, list(ind))
     }

     sq<-rep(1:len, coln*rown)

     if (all(rep.col.resp==0) && all(rep.row.resp==0)) return(NULL)

     for (i in 1:number.sample) {
          if (rep.row.resp[i]==1 && rep.col.resp[i]==1) ret[i]<-ifelse(ind.binary[i]==1,
                                                      rbinom(1, 1, sens.ind),
                                                      rbinom(1, 1, 1-spec.ind))
     }

     list(dframe=data.frame(x=cova, col.resp=rep.col.resp, row.resp=rep.row.resp,
          coln=col.groupn, rown=row.groupn, sqn=sq, retest=ret), ind=individual, prob=pij)

}

