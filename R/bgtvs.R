

"bgtvs" <- function
(n, s, Y, conf.level = 0.95, alternative="two.sided", maxiter=100)
{

alternative<-match.arg(alternative, choices=c("two.sided", "less", "greater"))

# error checking:

if (length(Y)!=length(n)||length(n)!=length(s))
 {stop("vectors n, s, Y must have exactly the same length")}

if ( any(Y>n) )
 {stop("values of Y must be smaller than or equal to the corresponding values n") }

if(conf.level<=0 | conf.level>=1)
 {stop("conf.level must be a numeric value between 0, and 1, usually e.g. 0.95")}



 m <- length(Y)
 alpha <- 1-conf.level

# # # Likelihood according to Hepworth 1996 (1):

likHepI <- function(n, s, Y, p)
{
 prod( choose(n,Y) * (( 1-(1-p)^s )^Y) * ( (1-p)^(s*(n-Y)) )  )
}

# function to estimate the MLE:

estBGTvgs2<-function(n, s, Y)
{
if (length(Y)!=length(n)||length(n)!=length(s)) {stop("vectors n, s, Y must have exactly the same length")}

# total number of individuals in the trial

# the current iteration method is sensitive 
# to give wrong results for too small tolerances

 tol=1e-10	
 total <- sum(n*s)

# Hepworth(1996), equation 2
 
mindiff <- function(n,s,Y,p, total)
{total - sum( (s*Y) / (1-(1-p)^(s)) )}

# check the case all Y=0 (equa 2 not defined),
# else iterate p

if(all(Y==0))
 {esti<-0}

 else
  {
   maxiter <- 100
   dir<-1
   crit<-numeric(length=maxiter)

    if(any(Y==0))
     {esti<-0+tol}
    else
     {esti<-0}

     crit[1] <- mindiff(n=n, s=s, Y=Y, p=esti, total=total)
     step=0.5

     for(i in 2:maxiter)
      {
       crit[i] <- mindiff(n=n, s=s, Y=Y, p=esti, total=total)
       if(sign(crit[i-1])!=sign(crit[i]))
        {dir <- dir*(-1); step <- step/2}
       if(esti>1 | esti<0)
        {dir<-dir*(-1)}
       esti<-esti+dir*step
      }
  }
return(esti)
}

# # # Point estimate for the observed outcome:
# Pobs (Hepworth)

point.esti <- estBGTvgs2(n=n, s=s, Y=Y)


# # # calculate the MLE for all possible combinations of Y1, ..., Ym

  possibleY<-list()
 
  for(i in 1:length(n))
  {
  possibleY[[i]]<-0:n[i]
  }
 allComb=expand.grid(possibleY)


 MLEComb<-numeric(length=nrow(allComb) )

 for (comb in 1:nrow(allComb))
  {Ytemp<-as.numeric(allComb[comb,])
   MLEComb[comb] <- estBGTvgs2(n=n, s=s, Y=Ytemp)
  }

 out=cbind(allComb,MLEComb)

# # # order the event space by their associated MLE

outlower<-as.matrix( out[which(out$MLEComb >= point.esti), 1:m] )

outupper <- as.matrix( out[which(out$MLEComb <= point.esti), 1:m])

# Iteration functions for the upper and lower bound:

itupperlimit <- function(pstart, outupper, alpha, n, s, maxiter=maxiter)
{
  dir<-1
  step<-0.1

   crit<-numeric(length=maxiter)

     crit[1] <- alpha - sum(apply(X=outupper, MARGIN=1, FUN=likHepI, n=n, s=s, p=pstart))
     
     p <- pstart+dir*step

     for(i in 2:maxiter)
      {
       crit[i] <- alpha - sum(apply(X=outupper, MARGIN=1, FUN=likHepI, n=n, s=s, p=p))
       if(sign(crit[i-1])!=sign(crit[i]))
        {dir <- dir*(-1); step <- step/2}

       p <- min(p+dir*step,1)
      }

return(pu=p)
}
   
 
itlowerlimit <- function(pstart, outlower, alpha, n, s, maxiter=maxiter)
{
  dir<-(-1)
  step<-0.1

   crit<-numeric(length=maxiter)

     crit[1] <- alpha - sum(apply(X=outlower, MARGIN=1, FUN=likHepI, n=n, s=s, p=pstart))
     
     p<-pstart+dir*step

     for(i in 2:maxiter)
      {
       crit[i] <- alpha - sum(apply(X=outlower, MARGIN=1, FUN=likHepI, n=n, s=s, p=p))
       if(sign(crit[i-1])!=sign(crit[i]))
        {dir <- dir*(-1); step <- step/2}

       p <- max(p+dir*step,0)
      }

return(pl=p)
}
   
# # # Calculation of Bounds

if(alternative=="two.sided")
 {
upper <- itupperlimit(pstart=point.esti, outupper=outupper, alpha=alpha/2, n=n, s=s, maxiter=100)
lower <- itlowerlimit(pstart=point.esti, outlower=outlower, alpha=alpha/2, n=n, s=s, maxiter=100)
 }

if(alternative=="less")
 {
 upper <- itupperlimit(pstart=point.esti, outupper=outupper, alpha=alpha, n=n, s=s, maxiter=100)
 lower <- 0
 }

if(alternative=="greater")
 {
 upper=1
 lower <- itlowerlimit(pstart=point.esti, outlower=outlower, alpha=alpha/2, n=n, s=s, maxiter=100)
 }

out <- list(conf.int=c(lower, upper),
estimate=point.esti, 
conf.level=conf.level,
alternative=alternative,
input=rbind("number of groups"=n,
"group size"=s,
"number of positive groups"=Y)
)

class(out)<-"bgtvs"

return(out)

}



"print.bgtvs"<-function
(x,...)
{
args<-list(...)

conf<-round(100*x$conf.level,2)

cat("Input","\n")
print(x$input)
cat("\n")
cat("Exact",conf,"percent confidence interval","\n")
cat("\n")

if(is.null(args$digits))
args$digits <- 4
args$x<-x[["conf.int"]]
do.call("print", args)
cat("\n")
cat("Point estimate:","\n")
cat("\n")
args$x<-x[["estimate"]]
do.call("print", args)
}
