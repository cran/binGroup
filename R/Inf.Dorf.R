
# Start  inf.dorf.measures() function
###################################################################
#    Brianna Hitt - 4-18-17
#    Purpose: calculates the testing expenditure and accuracy measures for informative Dorfman testing,
#             given a testing configuration; Informative Dorfman testing is the same as pool-specific
#             optimal Dorfman testing, described by McMahan et al. (2012a), except for that it attempts
#             to find the optimal testing configuration by examining all possible configurations rather
#             than using the greedy algorithm proposed by McMahan et al. (2012a)
#      calls: characteristics.pool() - calculates the testing expenditure for a given configuration
#             accuracy.dorf() - calculates measures of testing accuracy for informative Dorfman testing
#      inputs: p = probability that an individual is infected, can be an overall probability
#                  of disease or a vector of individual probabilities
#              Se = sensitivity of the diagnostic test
#              Sp = specificity of the diagnostic test
#              N = block size/initial group  size, up to 50
#              pool.sizes = a configuration/set of pool sizes from a matrix of all
#                    possible configurations, generated using the three-stage hierarchical setup
#      outputs: list of the testing expenditure (e), testing variance (v), and
#               measures of testing accuracy (summary)
#      notes: much of the code for this function, including the characteristics.pool() and
#             accuracy.dorf() functions are from Chris McMahan's programs, provided with
#             "Informative Dorfman screening" by McMahan, Tebbs, and Bilder (2012a)

inf.dorf.measures <- function(prob, se, sp, N, pool.sizes){

  # Saves original ordering
  ind.order<-order(prob)

  # Orders subjects, required under all Informative algorithms
  prob<-sort(prob)

  # Determines the number of subjects being screened and sets up vectors for storing summary measures
  N<-length(prob)
  pool.id<-rep(-1000,N)
  PSe<-rep(-100,N)
  PSp<-rep(-100,N)
  PPV<-rep(-100,N)
  NPV<-rep(-100,N)

  # pool.sizes is a single configuration/set of pool sizes from the matrix of all
  # possible configurations from the three-stage hierarchical testing setup
  psz <- pool.sizes
  J <- length(psz)

  # Finds measures pool by pool
  psz<-c(psz,0)
  lower<-1
  upper<-psz[1]
  vec.e<-rep(-1000,J)
  vec.v<-rep(-1000,J)
  for(i in 1:J){
    p.pool<-prob[lower:upper]
    pool.id[lower:upper]<-rep(i,length(p.pool))

    # calculates the testing expenditure and variance per individual
    res<-characteristics.pool(p=p.pool,se=se,sp=sp)
    vec.e[i]<-res$e
    vec.v[i]<-res$v

    # calculates measures of testing accuracy per individual
    res.acc<-accuracy.dorf(p=p.pool,se=se,sp=sp)
    PSe[lower:upper]<-res.acc$PSe
    PSp[lower:upper]<-res.acc$PSp
    PPV[lower:upper]<-res.acc$PPV
    NPV[lower:upper]<-res.acc$NPV

    lower<-1+upper
    upper<-upper+psz[i+1]
  }

  # Finds total expectation and variation
  res.e<-sum(vec.e)
  res.v<-sum(vec.v)

  # Returns all subjects to original ordering, along with their corresponding measures
  prob<-prob[order(ind.order)]
  pool.id<-pool.id[order(ind.order)]
  PSe<-PSe[order(ind.order)]
  PSp<-PSp[order(ind.order)]
  PPV<-PPV[order(ind.order)]
  NPV<-NPV[order(ind.order)]

  res.mat<-matrix(c(pool.id, prob, PSe, PSp, PPV, NPV),nrow=N,ncol=6,byrow=FALSE,
                  dimnames=list(as.character(1:N) , c("pool","probability","PSe", "PSp", "PPV", "NPV")))

  prob<-prob[ind.order]

  list("e"=res.e, "v"=res.v, "summary"=res.mat)
}

###################################################################




# Start  Inf.Dorf() function
###################################################################

#' @title Informative two-stage hierarchical (Dorfman) testing
#'
#' @author Brianna Hitt
#'
#' @description Find the optimal testing configuration for informative
#' two-stage hierarchical (Dorfman) testing and calculate the associated
#' operating characteristics.
#'
#' @details This function finds the optimal testing configuration and computes
#' the associated operating characteristics for informative two-stage
#' hierarchical (Dorfman) testing, implemented via the pool-specific optimal
#' Dorfman (PSOD) method described in McMahan et al. (2012a). This function
#' finds the optimal testing configuration by considering all possible configurations
#' instead of using the greedy algorithm proposed for PSOD testing.
#' See Hitt et al. (2018) at \url{http://www.chrisbilder.com/grouptesting/HBTM/}
#' McMahan et al. (2012a), or Dorfman (1943) for additional details.
#'
#' @param p the probability of disease, which can be an overall probability of disease,
#' from which a heterogeneous vector of individual probabilities will be generated, or
#' a vector of individual probabilities specified by the user
#' @param Se the sensitivity of the diagnostic test
#' @param Sp the specificity of the diagnostic test
#' @param group.sz a single block size for which to find the optimal configuration
#' out of all possible configurations and calculate the associated operating
#' characteristics, or a range of block sizes over which to find the optimal
#' testing configuration
#' @param obj.fn a list of objective functions which are minimized to find the
#' optimal testing configuration. Available options include "\kbd{ET}" (the expected
#' number of tests per individual), "\kbd{MAR}" (the expected number of tests divided
#' by the expected number of correct classifications, described in Malinovsky et al. (2016)),
#’ and "\kbd{GR}" (a linear combination of the expected number of tests, the number of
#’ misclassified negatives, and the number of misclassified positives, described in Graff &
#’ Roeloffs (1972)). See Hitt et al. (2018) at
#' \url{http://www.chrisbilder.com/grouptesting/HBTM/} for additional details.
#' @param weights a matrix of up to six sets of weights for the GR function. Each set of
#' weights is specified by a row of the matrix.
#' @param alpha a scale parameter for the beta distribution that specifies the degree of
#' heterogeneity for the generated probability vector
#'
#' @return A list containing:
#' \item{prob}{the overall probability of disease or vector of individual probabilities,
#' as specified by the user}
#' \item{alpha}{the level of heterogeneity used to generate the vector of individual
#' probabilities}
#' \item{Se}{the sensitivity of the diagnostic test}
#' \item{Sp}{the specificity of the diagnostic test}
#' \item{opt.ET, opt.MAR, opt.GR}{a list for each objective function specified by the
#' user, containing:
#' \describe{
#' \item{OTC}{a list specifying the optimal testing configuration, which includes:
#' \describe{
#' \item{Block.sz}{the block size/overall group size, which is not tested}
#' \item{pool.szs}{pool sizes for the first stage of testing}}}
#' \item{p.vec}{the sorted vector of individual probabilities}
#' \item{ET}{the expected testing expenditure for the OTC}
#' \item{value}{the value of the objective function per individual}
#' \item{PSe}{the overall pooling sensitivity for the algorithm}
#' \item{PSp}{the overall pooling specificity for the algorithm}
#' \item{PPPV}{the overall pooling positive predictive value for the algorithm}
#' \item{PNPV}{the overall pooling negative predictive value for the algorithm}}}
#'
#' @examples
#' # Find the optimal testing configuration for informative two-stage
#' #   hierarchical (Dorfman) testing
#' Inf.Dorf(p=0.01, Se=0.95, Sp=0.95, group.sz=3:20, obj.fn=c("ET", "MAR"),
#' weights=NULL, alpha=2)
#'
#' # Find the optimal testing configuration for informative two-stage
#' #  hierarchical (Dorfman) testing, for a specified vector of
#' #  individual probabilities
#' set.seed(8791)
#' Inf.Dorf(p=rbeta(10,2,200), Se=0.90, Sp=0.90, group.sz=10,
#' obj.fn=c("ET", "MAR", "GR"),
#' weights=matrix(data=c(1,1,10,10), nrow=2, ncol=2, byrow=TRUE),
#' alpha=NA)
#'
#' @seealso
#' \code{\link{NI.Dorf}} for non-informative
#' two-stage hierarchical (Dorfman) testing
#'
#' \code{\link{OTC}} for finding the
#' optimal testing configuration for a number of standard group testing algorithms
#'
#' \url{http://chrisbilder.com/grouptesting/HBTM}
#'
#' @references
#' \emph{Dorfman, R. (1943)}. The detection of defective members of large populations.
#' \emph{The Annals of Mathematical Statistics, 14, 436-440}.
#'
#' \emph{Graff, L.E. & Roeloffs, R. (1972)}. Group testing in the presence of test error;
#' an extension of the Dorfman procedure. \emph{Technometrics, 14, 113-122}.
#'
#' \emph{Malinovsky, Y.; Albert, P.S. & Roy, A. (2016)}. Reader reaction: A note on the
#' evaluation of group testing algorithms in the presence of misclassification.
#' \emph{Biometrics, 72, 299-302}.
#'
#' \emph{McMahan, C.S.; Tebbs, J.M. & Bilder, C.R. (2012a)}. Informative Dorfman screening.
#' \emph{Biometrics, 68, 287-296}.
#'
#' @family Optimal Testing Configuration functions

#    Brianna Hitt - 4-18-17
Inf.Dorf <- function(p, Se, Sp, group.sz, obj.fn, weights, alpha){

  start.time<-base::proc.time()

  set.of.blocks <- group.sz

  save.ET <- matrix(data = NA, nrow = length(set.of.blocks), ncol = max(set.of.blocks)+length(p)+10)
  save.MAR <- matrix(data = NA, nrow = length(set.of.blocks), ncol = max(set.of.blocks)+length(p)+10)
  save.GR1 <- matrix(data = NA, nrow = length(set.of.blocks), ncol = max(set.of.blocks)+length(p)+10)
  save.GR2 <- matrix(data = NA, nrow = length(set.of.blocks), ncol = max(set.of.blocks)+length(p)+10)
  save.GR3 <- matrix(data = NA, nrow = length(set.of.blocks), ncol = max(set.of.blocks)+length(p)+10)
  save.GR4 <- matrix(data = NA, nrow = length(set.of.blocks), ncol = max(set.of.blocks)+length(p)+10)
  save.GR5 <- matrix(data = NA, nrow = length(set.of.blocks), ncol = max(set.of.blocks)+length(p)+10)
  save.GR6 <- matrix(data = NA, nrow = length(set.of.blocks), ncol = max(set.of.blocks)+length(p)+10)

  count <- 1

  for(N in set.of.blocks){
    # build a vector of probabilities for a heterogeneous population
    if(length(p)==1){
      p.vec <- p.vec.func(p = p, alpha = alpha, grp.sz = N)
    } else if(length(p)>1){
      p.vec <- sort(p)
    }

    # generate a matrix of all possible configurations/sets of pool sizes
    # the partitions::parts() function requires loading the partitions library
    # do not include the first column, which would test the initial group twice
    possible.groups <- partitions::parts(n = N)[,-1]
    # could also limit the pool sizes to a maximum pool size of M
    #groups.M <- which(apply(possible.groups, 2, function(x) any(x>M)))
    #possible.groups <- possible.groups[,c(1,groups.M)]

    save.it <- matrix(data = NA, nrow = ncol(possible.groups), ncol = max(set.of.blocks)+length(p)+17)

    counter <- 1
    for(c in 1:ncol(possible.groups)){
      pool.sizes <- (possible.groups[,c])[possible.groups[,c]!=0]
      # calculate descriptive measures for informative Dorfman testing, given a configuration/set of pool sizes
      save.info <- inf.dorf.measures(prob = p.vec, se = Se, sp = Sp, N = N, pool.sizes = pool.sizes)

      # extract the configuration/pool sizes
      row.names(save.info$summary)=NULL
      pool.sz <- table(save.info$summary[,1])
      row.names(pool.sz)=NULL

      # extract accuracy measures for each individual
      ET <- save.info$e
      PSe.vec <- save.info$summary[,3]
      PSp.vec <- save.info$summary[,4]
      MAR <- MAR.func(ET, PSe.vec, PSp.vec, p.vec)

      # calculate overall accuracy measures
      PSe <- sum(p.vec*PSe.vec)/sum(p.vec)
      PSp <- sum((1-p.vec)*(PSp.vec))/sum(1-p.vec)
      PPPV <- sum(p.vec*PSe.vec)/sum(p.vec*PSe.vec + (1-p.vec)*(1-PSp.vec))
      PNPV <- sum((1-p.vec)*PSp.vec)/sum((1-p.vec)*PSp.vec + p.vec*(1-PSe.vec))

      # for each row in the matrix of weights, calculate the GR function
      if(is.null(dim(weights))){
        GR1 <- NA
        GR2 <- NA
        GR3 <- NA
        GR4 <- NA
        GR5 <- NA
        GR6 <- NA
      } else{
        GR1 <- GR.func(ET, p.vec, PSe.vec, PSp.vec, D1=weights[1,1], D2=weights[1,2])
        if(dim(weights)[1]>=2){
          GR2 <- GR.func(ET, p.vec, PSe.vec, PSp.vec, D1=weights[2,1], D2=weights[2,2])
        } else{GR2 <- NA}
        if(dim(weights)[1]>=3){
          GR3 <- GR.func(ET, p.vec, PSe.vec, PSp.vec, D1=weights[3,1], D2=weights[3,2])
        } else{GR3 <- NA}
        if(dim(weights)[1]>=4){
          GR4 <- GR.func(ET, p.vec, PSe.vec, PSp.vec, D1=weights[4,1], D2=weights[4,2])
        } else{GR4 <- NA}
        if(dim(weights)[1]>=5){
          GR5 <- GR.func(ET, p.vec, PSe.vec, PSp.vec, D1=weights[5,1], D2=weights[5,2])
        } else{GR5 <- NA}
        if(dim(weights)[1]>=6){
          GR6 <- GR.func(ET, p.vec, PSe.vec, PSp.vec, D1=weights[6,1], D2=weights[6,2])
        } else{GR6 <- NA}
      }

      save.it[counter,] <- c(sort(p), alpha, Se, Sp, N, ET, ET/N, MAR, GR1/N, GR2/N, GR3/N, GR4/N, GR5/N, GR6/N, PSe, PSp, PPPV, PNPV, pool.sz, rep(0, max(0, max(set.of.blocks)-length(pool.sz))))
      counter <- counter + 1

    }

    # find the best configuration for each block size N, out of all possible configurations
    save.ET[count,] <- save.it[save.it[,(length(p)+6)]==min(save.it[,(length(p)+6)]),c(1:(length(p)+5),(length(p)+6),(length(p)+14):ncol(save.it))]
    save.MAR[count,] <- save.it[save.it[,(length(p)+7)]==min(save.it[,(length(p)+7)]),c(1:(length(p)+5),(length(p)+7),(length(p)+14):ncol(save.it))]
    if(class(try(save.GR1[count,] <- save.it[save.it[,(length(p)+8)]==min(save.it[,(length(p)+8)]),c(1:(length(p)+5),(length(p)+8),(length(p)+14):ncol(save.it))],silent=T))!="try-error"){
      save.GR1[count,] <- save.it[save.it[,(length(p)+8)]==min(save.it[,(length(p)+8)]),c(1:(length(p)+5),(length(p)+8),(length(p)+14):ncol(save.it))]
    }
    if(class(try(save.GR2[count,] <- save.it[save.it[,(length(p)+9)]==min(save.it[,(length(p)+9)]),c(1:(length(p)+5),(length(p)+9),(length(p)+14):ncol(save.it))],silent=T))!="try-error"){
      save.GR2[count,] <- save.it[save.it[,(length(p)+9)]==min(save.it[,(length(p)+9)]),c(1:(length(p)+5),(length(p)+9),(length(p)+14):ncol(save.it))]
    }
    if(class(try(save.GR3[count,] <- save.it[save.it[,(length(p)+10)]==min(save.it[,(length(p)+10)]),c(1:(length(p)+5),(length(p)+10),(length(p)+14):ncol(save.it))],silent=T))!="try-error"){
      save.GR3[count,] <- save.it[save.it[,(length(p)+10)]==min(save.it[,(length(p)+10)]),c(1:(length(p)+5),(length(p)+10),(length(p)+14):ncol(save.it))]
    }
    if(class(try(save.GR4[count,] <- save.it[save.it[,(length(p)+11)]==min(save.it[,(length(p)+11)]),c(1:(length(p)+5),(length(p)+11),(length(p)+14):ncol(save.it))],silent=T))!="try-error"){
      save.GR4[count,] <- save.it[save.it[,(length(p)+11)]==min(save.it[,(length(p)+11)]),c(1:(length(p)+5),(length(p)+11),(length(p)+14):ncol(save.it))]
    }
    if(class(try( save.GR5[count,] <- save.it[save.it[,(length(p)+12)]==min(save.it[,(length(p)+12)]),c(1:(length(p)+5),(length(p)+12),(length(p)+14):ncol(save.it))],silent=T))!="try-error"){
      save.GR5[count,] <- save.it[save.it[,(length(p)+12)]==min(save.it[,(length(p)+12)]),c(1:(length(p)+5),(length(p)+12),(length(p)+14):ncol(save.it))]
    }
    if(class(try(save.GR6[count,] <- save.it[save.it[,(length(p)+13)]==min(save.it[,(length(p)+13)]),c(1:(length(p)+5),(length(p)+13),(length(p)+14):ncol(save.it))],silent=T))!="try-error"){
      save.GR6[count,] <- save.it[save.it[,(length(p)+13)]==min(save.it[,(length(p)+13)]),c(1:(length(p)+5),(length(p)+13),(length(p)+14):ncol(save.it))]
    }

    cat("Block Size =", N, "\n")
    count <- count + 1
  }

  # find the optimal configuration over all block sizes considered
  result.ET <- save.ET[save.ET[,(length(p)+6)]==min(save.ET[,(length(p)+6)]),]
  result.MAR <- save.MAR[save.MAR[,(length(p)+6)]==min(save.MAR[,(length(p)+6)]),]
  result.GR1 <- save.GR1[save.GR1[,(length(p)+6)]==min(save.GR1[,(length(p)+6)]),]
  result.GR2 <- save.GR2[save.GR2[,(length(p)+6)]==min(save.GR2[,(length(p)+6)]),]
  result.GR3 <- save.GR3[save.GR3[,(length(p)+6)]==min(save.GR3[,(length(p)+6)]),]
  result.GR4 <- save.GR4[save.GR4[,(length(p)+6)]==min(save.GR4[,(length(p)+6)]),]
  result.GR5 <- save.GR5[save.GR5[,(length(p)+6)]==min(save.GR5[,(length(p)+6)]),]
  result.GR6 <- save.GR6[save.GR6[,(length(p)+6)]==min(save.GR6[,(length(p)+6)]),]

  if(length(p)==1){
    p.vec.ET <- p.vec.func(p=result.ET[1], alpha=result.ET[2], grp.sz=result.ET[5])
    p.vec.MAR <- p.vec.func(p=result.MAR[1], alpha=result.MAR[2], grp.sz=result.MAR[5])
    p.vec.GR1 <- p.vec.func(p=result.GR1[1], alpha=result.GR1[2], grp.sz=result.GR1[5])
    p.vec.GR2 <- p.vec.func(p=result.GR2[1], alpha=result.GR2[2], grp.sz=result.GR2[5])
    p.vec.GR3 <- p.vec.func(p=result.GR3[1], alpha=result.GR3[2], grp.sz=result.GR3[5])
    p.vec.GR4 <- p.vec.func(p=result.GR4[1], alpha=result.GR4[2], grp.sz=result.GR4[5])
    p.vec.GR5 <- p.vec.func(p=result.GR5[1], alpha=result.GR5[2], grp.sz=result.GR5[5])
    p.vec.GR6 <- p.vec.func(p=result.GR6[1], alpha=result.GR6[2], grp.sz=result.GR6[5])
  } else if(length(p)>1){
    p.vec.ET <- result.ET[1:length(p)]
    p.vec.MAR <- result.MAR[1:length(p)]
    p.vec.GR1 <- result.GR1[1:length(p)]
    p.vec.GR2 <- result.GR2[1:length(p)]
    p.vec.GR3 <- result.GR3[1:length(p)]
    p.vec.GR4 <- result.GR4[1:length(p)]
    p.vec.GR5 <- result.GR5[1:length(p)]
    p.vec.GR6 <- result.GR6[1:length(p)]
  }

  # create a list of results for each objective function
  opt.ET <- list("OTC" = list("Block.sz" = result.ET[(length(p)+4)], "pool.szs" = (result.ET[(length(p)+11):length(result.ET)])[result.ET[(length(p)+11):length(result.ET)]!=0]), "p.vec" = p.vec.ET,
                 "ET" = result.ET[(length(p)+5)], "value" = result.ET[(length(p)+6)], "PSe" = result.ET[(length(p)+7)], "PSp" = result.ET[(length(p)+8)], "PPPV" = result.ET[(length(p)+9)], "PNPV" = result.ET[(length(p)+10)])
  opt.MAR <- list("OTC" = list("Block.sz" = result.MAR[(length(p)+4)], "pool.szs" = (result.MAR[(length(p)+11):length(result.MAR)])[result.MAR[(length(p)+11):length(result.MAR)]!=0]), "p.vec" = p.vec.MAR,
                  "ET" = result.MAR[(length(p)+5)], "value" = result.MAR[(length(p)+6)], "PSe" = result.MAR[(length(p)+7)], "PSp" = result.MAR[(length(p)+8)], "PPPV" = result.MAR[(length(p)+9)], "PNPV" = result.MAR[(length(p)+10)])
  opt.GR1 <- list("OTC" = list("Block.sz" = result.GR1[(length(p)+4)], "pool.szs" = (result.GR1[(length(p)+11):length(result.GR1)])[result.GR1[(length(p)+11):length(result.GR1)]!=0]), "p.vec" = p.vec.GR1,
                  "ET" = result.GR1[(length(p)+5)], "value" = result.GR1[(length(p)+6)], "PSe" = result.GR1[(length(p)+7)], "PSp" = result.GR1[(length(p)+8)], "PPPV" = result.GR1[(length(p)+9)], "PNPV" = result.GR1[(length(p)+10)])
  opt.GR2 <- list("OTC" = list("Block.sz" = result.GR2[(length(p)+4)], "pool.szs" = (result.GR2[(length(p)+11):length(result.GR2)])[result.GR2[(length(p)+11):length(result.GR2)]!=0]), "p.vec" = p.vec.GR2,
                  "ET" = result.GR2[(length(p)+5)], "value" = result.GR2[(length(p)+6)], "PSe" = result.GR2[(length(p)+7)], "PSp" = result.GR2[(length(p)+8)], "PPPV" = result.GR2[(length(p)+9)], "PNPV" = result.GR2[(length(p)+10)])
  opt.GR3 <- list("OTC" = list("Block.sz" = result.GR3[(length(p)+4)], "pool.szs" = (result.GR3[(length(p)+11):length(result.GR3)])[result.GR3[(length(p)+11):length(result.GR3)]!=0]), "p.vec" = p.vec.GR3,
                  "ET" = result.GR3[(length(p)+5)], "value" = result.GR3[(length(p)+6)], "PSe" = result.GR3[(length(p)+7)], "PSp" = result.GR3[(length(p)+8)], "PPPV" = result.GR3[(length(p)+9)], "PNPV" = result.GR3[(length(p)+10)])
  opt.GR4 <- list("OTC" = list("Block.sz" = result.GR4[(length(p)+4)], "pool.szs" = (result.GR4[(length(p)+11):length(result.GR4)])[result.GR4[(length(p)+11):length(result.GR4)]!=0]), "p.vec" = p.vec.GR4,
                  "ET" = result.GR4[(length(p)+5)], "value" = result.GR4[(length(p)+6)], "PSe" = result.GR4[(length(p)+7)], "PSp" = result.GR4[(length(p)+8)], "PPPV" = result.GR4[(length(p)+9)], "PNPV" = result.GR4[(length(p)+10)])
  opt.GR5 <- list("OTC" = list("Block.sz" = result.GR5[(length(p)+4)], "pool.szs" = (result.GR5[(length(p)+11):length(result.GR5)])[result.GR5[(length(p)+11):length(result.GR5)]!=0]), "p.vec" = p.vec.GR5,
                  "ET" = result.GR5[(length(p)+5)], "value" = result.GR5[(length(p)+6)], "PSe" = result.GR5[(length(p)+7)], "PSp" = result.GR5[(length(p)+8)], "PPPV" = result.GR5[(length(p)+9)], "PNPV" = result.GR5[(length(p)+10)])
  opt.GR6 <- list("OTC" = list("Block.sz" = result.GR6[(length(p)+4)], "pool.szs" = (result.GR6[(length(p)+11):length(result.GR6)])[result.GR6[(length(p)+11):length(result.GR6)]!=0]), "p.vec" = p.vec.GR6,
                  "ET" = result.GR6[(length(p)+5)], "value" = result.GR6[(length(p)+6)], "PSe" = result.GR6[(length(p)+7)], "PSp" = result.GR6[(length(p)+8)], "PPPV" = result.GR6[(length(p)+9)], "PNPV" = result.GR6[(length(p)+10)])

  # create a list of results, including all objective functions
  opt.all <- list("opt.ET" = opt.ET, "opt.MAR" = opt.MAR, "opt.GR1" = opt.GR1, "opt.GR2" = opt.GR2,
                  "opt.GR3" = opt.GR3, "opt.GR4" = opt.GR4, "opt.GR5" = opt.GR5, "opt.GR6" = opt.GR6)
  # remove any objective functions not requested by the user
  opt.req <- Filter(function(x) !is.na(x$ET), opt.all)
  time.it(start.time)
  c("prob" = list(p), "alpha" = alpha, "Se" = Se, "Sp" = Sp, opt.req)

}

###################################################################
