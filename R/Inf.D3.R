
# Start  Inf.D3() function
###################################################################

#' @title Informative three-stage hierarchical testing
#'
#' @author Brianna Hitt
#'
#' @description Find the optimal testing configuration for informative
#' three-stage hierarchical testing and calculate the associated
#' operating charcteristics.
#'
#' @details This function finds the optimal testing configuration and
#' computes the associated operating characteristics for informative
#' three-stage hierarchical testing. This function finds the optimal
#' testing configuration by considering all possible configurations.
#' See Hitt et al. (2018) at \url{http://www.chrisbilder.com/grouptesting/HBTM/}
#' or Black et al. (2015) for additional details.
#'
#' @param p the probability of disease, which can be an overall probability of disease,
#' from which a heterogeneous vector of individual probabilities will be generated, or
#' a vector of individual probabilities specified by the user
#' @param Se the sensitivity of the diagnostic test
#' @param Sp the specificity of the diagnostic test
#' @param group.sz a single group size over which to find the optimal testing
#' configuration out of all possible configurations, or a range of group sizes
#' over which to find the optimal testing configuration
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
#' \item{prob}{the overall probability of disease or vector of individual probabilities, as specified by the user}
#' \item{alpha}{the level of heterogeneity used to generate the vector of individual probabilities}
#' \item{Se}{the sensitivity of the diagnostic test}
#' \item{Sp}{the specificity of the diagnostic test}
#' \item{opt.ET, opt.MAR, opt.GR}{a list for each objective function specified by the user, containing:
#' \describe{
#' \item{OTC}{a list specifying the optimal testing configuration, which includes:
#' \describe{
#' \item{Stage1}{the pool size for the first stage of testing, i.e. the initial group size}
#' \item{Stage2}{pool sizes for the second stage of testing}}}
#' \item{p.vec}{the sorted vector of individual probabilities}
#' \item{ET}{the expected testing expenditure for the OTC}
#' \item{value}{the value of the objective function per individual}
#' \item{PSe}{the overall pooling sensitivity for the algorithm}
#' \item{PSp}{the overall pooling specificity for the algorithm}
#' \item{PPPV}{the overall pooling positive predictive value for the algorithm}
#' \item{PNPV}{the overall pooling negative predictive value for the algorithm}}}
#'
#' @examples
#' # Find the optimal testing configuration for information three-stage
#' #   hierarchical testing
#' Inf.D3(p=0.05, Se=0.99, Sp=0.99, group.sz=10:15, obj.fn=c("ET", "MAR"),
#' weights=NULL, alpha=0.5)
#'
#' # Find the optimal testing configuration out of all possible configurations
#' #   for a specified group size and vector of individual probabilities
#' set.seed(82763)
#' Inf.D3(p=rbeta(12,2,200), Se=0.99, Sp=0.99, group.sz=12,
#' obj.fn=c("ET", "MAR", "GR"),
#' weights=matrix(data=c(1,1), nrow=1, ncol=2, byrow=TRUE),
#' alpha=NA)
#'
#' @seealso
#' \code{\link{NI.D3}} for non-informative three-stage hierarchical testing
#'
#' \code{\link{OTC}} for finding the optimal testing configuration for a
#' number of standard group testing algorithms
#'
#' \url{http://chrisbilder.com/grouptesting/HBTM}
#'
#' @references
#' \emph{Black, M.S.; Bilder, C.R. & Tebbs, J.M. (2015)}.
#' Optimal retesting configurations for hierarchical group testing.
#' \emph{Journal of the Royal Statistical Society. Series C: Applied Statistics,
#' 64, 693-710}.
#'
#' \emph{Graff, L.E. & Roeloffs, R. (1972)}. Group testing in the presence of
#' test error; an extension of the Dorfman procedure. \emph{Technometrics, 14, 113-122}.
#'
#' \emph{Malinovsky, Y.; Albert, P.S. & Roy, A. (2016)}. Reader reaction: A note on the
#' evaluation of group testing algorithms in the presence of misclassification.
#' \emph{Biometrics, 72, 299-302}.
#'
#' @family Optimal Testing Configuration functions

#    Brianna Hitt - 05-01-17

Inf.D3 <- function(p, Se, Sp, group.sz, obj.fn, weights, alpha){

  start.time<-base::proc.time()

  set.of.I <- group.sz

  save.ET <- matrix(data = NA, nrow = length(set.of.I), ncol = max(set.of.I)+length(p)+10)
  save.MAR <- matrix(data = NA, nrow = length(set.of.I), ncol = max(set.of.I)+length(p)+10)
  save.GR1 <- matrix(data = NA, nrow = length(set.of.I), ncol = max(set.of.I)+length(p)+10)
  save.GR2 <- matrix(data = NA, nrow = length(set.of.I), ncol = max(set.of.I)+length(p)+10)
  save.GR3 <- matrix(data = NA, nrow = length(set.of.I), ncol = max(set.of.I)+length(p)+10)
  save.GR4 <- matrix(data = NA, nrow = length(set.of.I), ncol = max(set.of.I)+length(p)+10)
  save.GR5 <- matrix(data = NA, nrow = length(set.of.I), ncol = max(set.of.I)+length(p)+10)
  save.GR6 <- matrix(data = NA, nrow = length(set.of.I), ncol = max(set.of.I)+length(p)+10)

  count <- 1

  for(I in set.of.I){
    # build a vector of probabilities for a heterogeneous population
    if(length(p)==1){
      p.vec <- p.vec.func(p = p, alpha = alpha, grp.sz = I)
    } else if(length(p)>1){
      p.vec <- sort(p)
    }

    # generate a matrix of all possible configurations/sets of pool sizes
    # the partitions::parts() function requires loading the partitions library
    # do not include the first column, which would test the initial group twice
    possible.groups <- partitions::parts(n = I)[,-1]
    # could also ignore all columns that include pool sizes of 1
    #groups.1 <- which(apply(possible.groups, 2, function(x) any(x==1)))
    #possible.groups <- possible.groups[,c(1,groups.1)]

    save.it <- matrix(data = NA, nrow = ncol(possible.groups), ncol = max(set.of.I)+length(p)+17)

    counter <- 1
    for(c in 1:ncol(possible.groups)){
      # extract the configuration, ordering, and group sizes for each column
      config <- c(possible.groups[,c], 1:I)
      order.for.p <- config[(1+I):(2*I)]
      gp.sizes <- config[1:I]

      # call hierarchical.desc2() for the configuration
      save.info <- hierarchical.desc2(p = p.vec[order.for.p],se = Se, sp = Sp, I2 = gp.sizes[gp.sizes!=0], order.p=FALSE)

      # extract accuracy measures for each individual
      ET <- save.info$ET
      PSe.vec <- save.info$individual.testerror$pse.vec
      PSp.vec <- save.info$individual.testerror$psp.vec
      MAR <- MAR.func(ET, PSe.vec, PSp.vec, p.vec)

      # calculate overall accuracy measures
      group.testerror <- save.info$group.testerror
      names(group.testerror) <- NULL
      PSe <- group.testerror[1]
      PSp <- group.testerror[2]
      PPPV <- group.testerror[3]
      PNPV <- group.testerror[4]

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

      save.it[counter,] <- c(sort(p), alpha, Se, Sp, I, ET, ET/I, MAR, GR1/I, GR2/I, GR3/I, GR4/I, GR5/I, GR6/I, PSe, PSp, PPPV, PNPV, gp.sizes, rep(0, max(0, max(set.of.I)-length(gp.sizes))))
      counter <- counter + 1

    }

    # find the best configuration for each initial group size I, out of all possible configurations
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

    cat("Initial Group Size =", I, "\n")
    count <- count + 1
  }

  # find the optimal testing configuration, over all initial group sizes considered
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
  opt.ET <- list("OTC" = list("Stage1" = result.ET[(length(p)+4)], "Stage2" = (result.ET[(length(p)+11):length(result.ET)])[result.ET[(length(p)+11):length(result.ET)]!=0]), "p.vec" = p.vec.ET,
                 "ET" = result.ET[(length(p)+5)], "value" = result.ET[(length(p)+6)], "PSe" = result.ET[(length(p)+7)], "PSp" = result.ET[(length(p)+8)], "PPPV" = result.ET[(length(p)+9)], "PNPV" = result.ET[(length(p)+10)])
  opt.MAR <- list("OTC" = list("Stage1" = result.MAR[(length(p)+4)], "Stage2" = (result.MAR[(length(p)+11):length(result.MAR)])[result.MAR[(length(p)+11):length(result.MAR)]!=0]), "p.vec" = p.vec.MAR,
                  "ET" = result.MAR[(length(p)+5)], "value" = result.MAR[(length(p)+6)], "PSe" = result.MAR[(length(p)+7)], "PSp" = result.MAR[(length(p)+8)], "PPPV" = result.MAR[(length(p)+9)], "PNPV" = result.MAR[(length(p)+10)])
  opt.GR1 <- list("OTC" = list("Stage1" = result.GR1[(length(p)+4)], "Stage2" = (result.GR1[(length(p)+11):length(result.GR1)])[result.GR1[(length(p)+11):length(result.GR1)]!=0]), "p.vec" = p.vec.GR1,
                  "ET" = result.GR1[(length(p)+5)], "value" = result.GR1[(length(p)+6)], "PSe" = result.GR1[(length(p)+7)], "PSp" = result.GR1[(length(p)+8)], "PPPV" = result.GR1[(length(p)+9)], "PNPV" = result.GR1[(length(p)+10)])
  opt.GR2 <- list("OTC" = list("Stage1" = result.GR2[(length(p)+4)], "Stage2" = (result.GR2[(length(p)+11):length(result.GR2)])[result.GR2[(length(p)+11):length(result.GR2)]!=0]), "p.vec" = p.vec.GR2,
                  "ET" = result.GR2[(length(p)+5)], "value" = result.GR2[(length(p)+6)], "PSe" = result.GR2[(length(p)+7)], "PSp" = result.GR2[(length(p)+8)], "PPPV" = result.GR2[(length(p)+9)], "PNPV" = result.GR2[(length(p)+10)])
  opt.GR3 <- list("OTC" = list("Stage1" = result.GR3[(length(p)+4)], "Stage2" = (result.GR3[(length(p)+11):length(result.GR3)])[result.GR3[(length(p)+11):length(result.GR3)]!=0]), "p.vec" = p.vec.GR3,
                  "ET" = result.GR3[(length(p)+5)], "value" = result.GR3[(length(p)+6)], "PSe" = result.GR3[(length(p)+7)], "PSp" = result.GR3[(length(p)+8)], "PPPV" = result.GR3[(length(p)+9)], "PNPV" = result.GR3[(length(p)+10)])
  opt.GR4 <- list("OTC" = list("Stage1" = result.GR4[(length(p)+4)], "Stage2" = (result.GR4[(length(p)+11):length(result.GR4)])[result.GR4[(length(p)+11):length(result.GR4)]!=0]), "p.vec" = p.vec.GR4,
                  "ET" = result.GR4[(length(p)+5)], "value" = result.GR4[(length(p)+6)], "PSe" = result.GR4[(length(p)+7)], "PSp" = result.GR4[(length(p)+8)], "PPPV" = result.GR4[(length(p)+9)], "PNPV" = result.GR4[(length(p)+10)])
  opt.GR5 <- list("OTC" = list("Stage1" = result.GR5[(length(p)+4)], "Stage2" = (result.GR5[(length(p)+11):length(result.GR5)])[result.GR5[(length(p)+11):length(result.GR5)]!=0]), "p.vec" = p.vec.GR5,
                  "ET" = result.GR5[(length(p)+5)], "value" = result.GR5[(length(p)+6)], "PSe" = result.GR5[(length(p)+7)], "PSp" = result.GR5[(length(p)+8)], "PPPV" = result.GR5[(length(p)+9)], "PNPV" = result.GR5[(length(p)+10)])
  opt.GR6 <- list("OTC" = list("Stage1" = result.GR6[(length(p)+4)], "Stage2" = (result.GR6[(length(p)+11):length(result.GR6)])[result.GR6[(length(p)+11):length(result.GR6)]!=0]), "p.vec" = p.vec.GR6,
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

