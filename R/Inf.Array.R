
# Start  Inf.Array() function
###################################################################

#' @title Informative array testing without master pooling
#'
#' @author Brianna Hitt
#'
#' @description Find the optimal testing configuration for
#' informative array testing without master pooling and calculate
#' the associated operating characteristics.
#'
#' @details This function finds the optimal testing configuration and
#' computes the associated operating characteristics for informative
#' array testing without master pooling, implemented using the
#' gradient arrangement described in McMahan et al. (2012b). Array
#' testing with master pooling involves amalgamating specimens in rows
#' and columns for the first stage of testing. This function uses only
#' square arrays, which is the way array-based group testing is carried
#' out in most real-world applications. See Hitt et al. (2018)
#' at \url{http://www.chrisbilder.com/grouptesting/HBTM/}
#' McMahan et al. (2012b) for additional details.
#'
#' @param p the probability of disease, which can be an overall probability of disease,
#' from which a heterogeneous vector of individual probabilities will be generated, or
#' a vector of individual probabilities specified by the user
#' @param Se the sensitivity of the diagnostic test
#' @param Sp the specificity of the diagnostic test
#' @param group.sz a single group size (representing the row/column size)
#' for which to calculate the operating characteristics, or a range of group
#' (row/column) sizes over which to find the optimal testing configuration
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
#' \item{Array.dim}{the row/column size for the first stage of testing}
#' \item{Array.sz}{the overall array size (the square of the row/column size)}}}
#' \item{p.mat}{the sorted matrix of individual probabilities, arranged using the gradient method described by McMahan et al. (2012b)}
#' \item{ET}{the expected testing expenditure for the OTC}
#' \item{value}{the value of the objective function per individual}
#' \item{PSe}{the overall pooling sensitivity for the algorithm}
#' \item{PSp}{the overall pooling specificity for the algorithm}
#' \item{PPPV}{the overall pooling positive predictive value for the algorithm}
#' \item{PNPV}{the overall pooling negative predictive value for the algorithm}}}
#'
#' @examples
#' # Find the optimal testing configuration for informative array
#' #   testing without master pooling over a range of group
#' #   (row/column) sizes
#' \dontrun{Inf.Array(p=0.03, Se=0.99, Sp=0.99, group.sz=5:10, obj.fn=c("ET", "MAR"),
#' weights=NULL, alpha=2)}
#'
#' # Find the optimal testing configuration for a specified group
#' #   (row/column) size for informative array testing without
#' #   master pooling
#' Inf.Array(p=rbeta(10,2,200), Se=0.95, Sp=0.95, group.sz=10,
#' obj.fn=c("ET", "MAR", "GR"),
#' weights=matrix(data=c(1,1,10,10), nrow=2, ncol=2, byrow=TRUE),
#' alpha=NA)
#'
#' @seealso
#' \code{\link{NI.Array}} for non-informative array testing
#' without master pooling
#'
#' \code{\link{NI.A2M}} for non-informative array testing
#' with master pooling
#'
#' \code{\link{OTC}} for finding the optimal testing configuration for
#' a number of standard group testing algorithms
#'
#' \url{http://chrisbilder.com/grouptesting/HBTM}
#'
#' @references
#' \emph{Graff, L.E. & Roeloffs, R. (1972)}. Group testing in the presence of
#' test error; an extension of the Dorfman procedure. \emph{Technometrics, 14, 113-122}.
#'
#' \emph{Malinovsky, Y.; Albert, P.S. & Roy, A. (2016)}. Reader reaction: A note on
#' the evaluation of group testing algorithms in the presence of misclassification.
#' \emph{Biometrics, 72, 299-302}.
#'
#' \emph{McMahan, C.S.; Tebbs, J.M. & Bilder, C.R. (2012b)}. Two-dimensional informative
#' array testing. \emph{Biometrics, 68, 793-804}.
#'
#' @family Optimal Testing Configuration functions

#    Brianna Hitt - 05-01-17

Inf.Array <- function(p, Se, Sp, group.sz, obj.fn, weights, alpha){

  start.time<-base::proc.time()

  set.of.I <- group.sz

  save.it <- matrix(data=NA, nrow=length(set.of.I), ncol=length(p)+18)
  count <- 1

  for(I in set.of.I){
    N <- I^2

    # build a vector of probabilities for a heterogeneous population
    if(length(p)==1){
      p.vec <- p.vec.func(p = p, alpha = alpha, grp.sz = N)
    } else if(length(p)>1){
      p.vec <- sort(p)
    }

    # build a matrix of probabilities using the gradient design
    p.ga <- Informative.array.prob(prob.vec = p.vec, nr = I, nc = I, method = "gd")

    # call Array.Measures() to calculate descriptive measures for the given array size
    save.info <- Array.Measures(p = p.ga, se = Se, sp = Sp)

    # extract accuracy measures for each individual
    ET <- save.info$T
    PSe.mat <- save.info$PSe
    PSp.mat <- save.info$PSp
    MAR <- MAR.func(ET, PSe.mat, PSp.mat, p.ga)

    # calculate overall accuracy measures
    PSe <- sum(p.ga*PSe.mat)/sum(p.ga)
    PSp <- sum((1-p.ga)*(PSp.mat))/sum(1-p.ga)
    PPPV <- sum(p.ga*PSe.mat)/sum(p.ga*PSe.mat + (1-p.ga)*(1-PSp.mat))
    PNPV <- sum((1-p.ga)*PSp.mat)/sum((1-p.ga)*PSp.mat + p.ga*(1-PSe.mat))

    # for each row in the matrix of weights, calculate the GR function
    if(is.null(dim(weights))){
      GR1 <- NA
      GR2 <- NA
      GR3 <- NA
      GR4 <- NA
      GR5 <- NA
      GR6 <- NA
    } else{
      GR1 <- GR.func(ET, p.ga, PSe.mat, PSp.mat, D1=weights[1,1], D2=weights[1,2])
      if(dim(weights)[1]>=2){
        GR2 <- GR.func(ET, p.ga, PSe.mat, PSp.mat, D1=weights[2,1], D2=weights[2,2])
      } else{GR2 <- NA}
      if(dim(weights)[1]>=3){
        GR3 <- GR.func(ET, p.ga, PSe.mat, PSp.mat, D1=weights[3,1], D2=weights[3,2])
      } else{GR3 <- NA}
      if(dim(weights)[1]>=4){
        GR4 <- GR.func(ET, p.ga, PSe.mat, PSp.mat, D1=weights[4,1], D2=weights[4,2])
      } else{GR4 <- NA}
      if(dim(weights)[1]>=5){
        GR5 <- GR.func(ET, p.ga, PSe.mat, PSp.mat, D1=weights[5,1], D2=weights[5,2])
      } else{GR5 <- NA}
      if(dim(weights)[1]>=6){
        GR6 <- GR.func(ET, p.ga, PSe.mat, PSp.mat, D1=weights[6,1], D2=weights[6,2])
      } else{GR6 <- NA}
    }

    save.it[count,] <- c(p, alpha, Se, Sp, I, N, ET, ET/N, MAR, GR1/N, GR2/N, GR3/N, GR4/N, GR5/N, GR6/N, PSe, PSp, PPPV, PNPV)
    cat("Row/Column Size = ", I, ", Array Size = ", N, "\n", sep="")
    count <- count + 1

  }

  # find the optimal testing configuration, over all array sizes considered
  result.ET <- save.it[save.it[,(length(p)+7)]==min(save.it[,(length(p)+7)]),c(1:(length(p)+6),(length(p)+7),(length(p)+15):ncol(save.it))]
  result.MAR <- save.it[save.it[,(length(p)+8)]==min(save.it[,(length(p)+8)]),c(1:(length(p)+6),(length(p)+8),(length(p)+15):ncol(save.it))]
  result.GR1 <- save.it[save.it[,(length(p)+9)]==min(save.it[,(length(p)+9)]),c(1:(length(p)+6),(length(p)+9),(length(p)+15):ncol(save.it))]
  result.GR2 <- save.it[save.it[,(length(p)+10)]==min(save.it[,(length(p)+10)]),c(1:(length(p)+6),(length(p)+10),(length(p)+15):ncol(save.it))]
  result.GR3 <- save.it[save.it[,(length(p)+11)]==min(save.it[,(length(p)+11)]),c(1:(length(p)+6),(length(p)+11),(length(p)+15):ncol(save.it))]
  result.GR4 <- save.it[save.it[,(length(p)+12)]==min(save.it[,(length(p)+12)]),c(1:(length(p)+6),(length(p)+12),(length(p)+15):ncol(save.it))]
  result.GR5 <- save.it[save.it[,(length(p)+13)]==min(save.it[,(length(p)+13)]),c(1:(length(p)+6),(length(p)+13),(length(p)+15):ncol(save.it))]
  result.GR6 <- save.it[save.it[,(length(p)+14)]==min(save.it[,(length(p)+14)]),c(1:(length(p)+6),(length(p)+14),(length(p)+15):ncol(save.it))]

  if(length(p)==1){
    p.vec.ET <- p.vec.func(p=result.ET[1], alpha=result.ET[2], grp.sz=result.ET[6])
    p.ga.ET <- Informative.array.prob(prob.vec = p.vec.ET, nr = result.ET[5], nc = result.ET[5], method = "gd")
    p.vec.MAR <- p.vec.func(p=result.MAR[1], alpha=result.MAR[2], grp.sz=result.MAR[6])
    p.ga.MAR <- Informative.array.prob(prob.vec = p.vec.MAR, nr = result.MAR[5], nc = result.MAR[5], method = "gd")
    if(is.null(dim(weights))){
      p.ga.GR1 <- NA
      p.ga.GR2 <- NA
      p.ga.GR3 <- NA
      p.ga.GR4 <- NA
      p.ga.GR5 <- NA
      p.ga.GR6 <- NA
    } else{
      p.vec.GR1 <- p.vec.func(p=result.GR1[1], alpha=result.GR1[2], grp.sz=result.GR1[6])
      p.ga.GR1 <- Informative.array.prob(prob.vec = p.vec.GR1, nr = result.GR1[5], nc = result.GR1[5], method = "gd")
      if(dim(weights)[1]>=2){
        p.vec.GR2 <- p.vec.func(p=result.GR2[1], alpha=result.GR2[2], grp.sz=result.GR2[6])
        p.ga.GR2 <- Informative.array.prob(prob.vec = p.vec.GR2, nr = result.GR2[5], nc = result.GR2[5], method = "gd")
      } else{p.ga.GR2 <- NA}
      if(dim(weights)[1]>=3){
        p.vec.GR3 <- p.vec.func(p=result.GR3[1], alpha=result.GR3[2], grp.sz=result.GR3[6])
        p.ga.GR3 <- Informative.array.prob(prob.vec = p.vec.GR3, nr = result.GR3[5], nc = result.GR3[5], method = "gd")
      } else{p.ga.GR3 <- NA}
      if(dim(weights)[1]>=4){
        p.vec.GR4 <- p.vec.func(p=result.GR4[1], alpha=result.GR4[2], grp.sz=result.GR4[6])
        p.ga.GR4 <- Informative.array.prob(prob.vec = p.vec.GR4, nr = result.GR4[5], nc = result.GR4[5], method = "gd")
      } else{p.ga.GR4 <- NA}
      if(dim(weights)[1]>=5){
        p.vec.GR5 <- p.vec.func(p=result.GR5[1], alpha=result.GR5[2], grp.sz=result.GR5[6])
        p.ga.GR5 <- Informative.array.prob(prob.vec = p.vec.GR5, nr = result.GR5[5], nc = result.GR5[5], method = "gd")
      } else{p.ga.GR5 <- NA}
      if(dim(weights)[1]>=6){
        p.vec.GR6 <- p.vec.func(p=result.GR6[1], alpha=result.GR6[2], grp.sz=result.GR6[6])
        p.ga.GR6 <- Informative.array.prob(prob.vec = p.vec.GR6, nr = result.GR6[5], nc = result.GR6[5], method = "gd")
      } else{p.ga.GR6 <- NA}
    }
  } else if(length(p)>1){
    p.vec <- sort(p)
    p.ga.ET <- Informative.array.prob(prob.vec = p.vec, nr = result.ET[(length(p)+4)], nc = result.ET[(length(p)+4)], method = "gd")
    p.ga.MAR <- Informative.array.prob(prob.vec = p.vec, nr = result.MAR[(length(p)+4)], nc = result.MAR[(length(p)+4)], method = "gd")
    if(is.null(dim(weights))){
      p.ga.GR1 <- NA
      p.ga.GR2 <- NA
      p.ga.GR3 <- NA
      p.ga.GR4 <- NA
      p.ga.GR5 <- NA
      p.ga.GR6 <- NA
    } else{
      p.ga.GR1 <- Informative.array.prob(prob.vec = p.vec, nr = result.GR1[(length(p)+4)], nc = result.GR1[(length(p)+4)], method = "gd")
      if(dim(weights)[1]>=2){
        p.ga.GR2 <- Informative.array.prob(prob.vec = p.vec, nr = result.GR2[(length(p)+4)], nc = result.GR2[(length(p)+4)], method = "gd")
      } else{p.ga.GR2 <- NA}
      if(dim(weights)[1]>=3){
        p.ga.GR3 <- Informative.array.prob(prob.vec = p.vec, nr = result.GR3[(length(p)+4)], nc = result.GR3[(length(p)+4)], method = "gd")
      } else{p.ga.GR3 <- NA}
      if(dim(weights)[1]>=4){
        p.ga.GR4 <- Informative.array.prob(prob.vec = p.vec, nr = result.GR4[(length(p)+4)], nc = result.GR4[(length(p)+4)], method = "gd")
      } else{p.ga.GR4 <- NA}
      if(dim(weights)[1]>=5){
        p.ga.GR5 <- Informative.array.prob(prob.vec = p.vec, nr = result.GR5[(length(p)+4)], nc = result.GR5[(length(p)+4)], method = "gd")
      } else{p.ga.GR5 <- NA}
      if(dim(weights)[1]>=6){
        p.ga.GR6 <- Informative.array.prob(prob.vec = p.vec, nr = result.GR6[(length(p)+4)], nc = result.GR6[(length(p)+4)], method = "gd")
      } else{p.ga.GR6 <- NA}
    }
  }

  # create a list of results for each objective function
  opt.ET <- list("OTC" = list("Array.dim" = result.ET[(length(p)+4)], "Array.sz" = result.ET[(length(p)+5)]), "p.mat" = p.ga.ET,
                 "ET" = result.ET[(length(p)+6)], "value" = result.ET[(length(p)+7)], "PSe" = result.ET[(length(p)+8)], "PSp" = result.ET[(length(p)+9)], "PPPV" = result.ET[(length(p)+10)], "PNPV" = result.ET[(length(p)+11)])
  opt.MAR <- list("OTC" = list("Array.dim" = result.MAR[(length(p)+4)], "Array.sz" = result.MAR[(length(p)+5)]), "p.mat" = p.ga.MAR,
                  "ET" = result.MAR[(length(p)+6)], "value" = result.MAR[(length(p)+7)], "PSe" = result.MAR[(length(p)+8)], "PSp" = result.MAR[(length(p)+9)], "PPPV" = result.MAR[(length(p)+10)], "PNPV" = result.MAR[(length(p)+11)])
  opt.GR1 <- list("OTC" = list("Array.dim" = result.GR1[(length(p)+4)], "Array.sz" = result.GR1[(length(p)+5)]), "p.mat" = p.ga.GR1,
                  "ET" = result.GR1[(length(p)+6)], "value" = result.GR1[(length(p)+7)], "PSe" = result.GR1[(length(p)+8)], "PSp" = result.GR1[(length(p)+9)], "PPPV" = result.GR1[(length(p)+10)], "PNPV" = result.GR1[(length(p)+11)])
  opt.GR2 <- list("OTC" = list("Array.dim" = result.GR2[(length(p)+4)], "Array.sz" = result.GR2[(length(p)+5)]), "p.mat" = p.ga.GR2,
                  "ET" = result.GR2[(length(p)+6)], "value" = result.GR2[(length(p)+7)], "PSe" = result.GR2[(length(p)+8)], "PSp" = result.GR2[(length(p)+9)], "PPPV" = result.GR2[(length(p)+10)], "PNPV" = result.GR2[(length(p)+11)])
  opt.GR3 <- list("OTC" = list("Array.dim" = result.GR3[(length(p)+4)], "Array.sz" = result.GR3[(length(p)+5)]), "p.mat" = p.ga.GR3,
                  "ET" = result.GR3[(length(p)+6)], "value" = result.GR3[(length(p)+7)], "PSe" = result.GR3[(length(p)+8)], "PSp" = result.GR3[(length(p)+9)], "PPPV" = result.GR3[(length(p)+10)], "PNPV" = result.GR3[(length(p)+11)])
  opt.GR4 <- list("OTC" = list("Array.dim" = result.GR4[(length(p)+4)], "Array.sz" = result.GR4[(length(p)+5)]), "p.mat" = p.ga.GR4,
                  "ET" = result.GR4[(length(p)+6)], "value" = result.GR4[(length(p)+7)], "PSe" = result.GR4[(length(p)+8)], "PSp" = result.GR4[(length(p)+9)], "PPPV" = result.GR4[(length(p)+10)], "PNPV" = result.GR4[(length(p)+11)])
  opt.GR5 <- list("OTC" = list("Array.dim" = result.GR5[(length(p)+4)], "Array.sz" = result.GR5[(length(p)+5)]), "p.mat" = p.ga.GR5,
                  "ET" = result.GR5[(length(p)+6)], "value" = result.GR5[(length(p)+7)], "PSe" = result.GR5[(length(p)+8)], "PSp" = result.GR5[(length(p)+9)], "PPPV" = result.GR5[(length(p)+10)], "PNPV" = result.GR5[(length(p)+11)])
  opt.GR6 <- list("OTC" = list("Array.dim" = result.GR6[(length(p)+4)], "Array.sz" = result.GR6[(length(p)+5)]), "p.mat" = p.ga.GR6,
                  "ET" = result.GR6[(length(p)+6)], "value" = result.GR6[(length(p)+7)], "PSe" = result.GR6[(length(p)+8)], "PSp" = result.GR6[(length(p)+9)], "PPPV" = result.GR6[(length(p)+10)], "PNPV" = result.GR6[(length(p)+11)])

  # create a list of results, including all objective functions
  opt.all <- list("opt.ET" = opt.ET, "opt.MAR" = opt.MAR, "opt.GR1" = opt.GR1, "opt.GR2" = opt.GR2,
                  "opt.GR3" = opt.GR3, "opt.GR4" = opt.GR4, "opt.GR5" = opt.GR5, "opt.GR6" = opt.GR6)
  # remove any objective functions not requested by the user
  opt.req <- Filter(function(x) !is.na(x$ET), opt.all)
  time.it(start.time)
  c("prob" = list(p), "alpha" = alpha, "Se" = Se, "Sp" = Sp, opt.req)

}

###################################################################
