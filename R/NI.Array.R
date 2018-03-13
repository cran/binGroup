
# Start  NI.Array() function
###################################################################

#' @title Non-informative array testing without master pooling
#'
#' @author Brianna Hitt
#'
#' @description Find the optimal testing configuration for
#' non-informative array testing without master pooling and
#' calculate the associated operating characteristics.
#'
#' @details This function finds the optimal testing configuration and
#' computes the associated operating characteristics for non-informative
#' array testing without master pooling. Array testing without master
#' pooling involves amalgamating specimens in rows and columns for the
#' first stage of testing. This function uses only square arrays,
#' which is the way array-based group testing is carried
#' out in most real-world applications. See Hitt et al. (2018) at
#' \url{http://www.chrisbilder.com/grouptesting/HBTM/} or
#' Kim et al. (2007) for additional details.
#'
#' @param p the probability of disease, which can be specified as an overall
#' probability of disease or a homogeneous vector of individual probabilities
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
#'
#' @return A list containing:
#' \item{prob}{the overall probability of disease or vector of individual probabilities,
#' as specified by the user}
#' \item{Se}{the sensitivity of the diagnostic test}
#' \item{Sp}{the specificity of the diagnostic test}
#' \item{opt.ET, opt.MAR, opt.GR}{a list for each objective function specified by the
#' user, containing:
#' \describe{
#' \item{OTC}{a list specifying the optimal testing configuration, which includes:
#' \describe{
#' \item{Array.dim}{the row/column size for the first stage of testing}
#' \item{Array.sz}{the overall array size (the square of the row/column size)}}}
#' \item{p.mat}{the matrix of individual probabilities}
#' \item{ET}{the expected testing expenditure for the OTC}
#' \item{value}{the value of the objective function per individual}
#' \item{PSe}{the overall pooling sensitivity for the algorithm}
#' \item{PSp}{the overall pooling specificity for the algorithm}
#' \item{PPPV}{the overall pooling positive predictive value for the algorithm}
#' \item{PNPV}{the overall pooling negative predictive value for the algorithm}}}
#'
#' @examples
#' # Find the optimal testing configuration for non-informative array
#' #   testing without master pooling over a range of group
#' #   (row/column) sizes
#' NI.Array(p=0.04, Se=0.95, Sp=0.95, group.sz=3:10,
#' obj.fn=c("ET", "MAR"), weights=NULL)
#'
#' # Find the optimal testing configuration for a specified group
#' #   (row/column) size for non-informative array testing without
#' #   master pooling
#' NI.Array(p=rep(0.01,8), Se=0.90, Sp=0.90, group.sz=8,
#' obj.fn=c("ET", "MAR", "GR"),
#' weights=matrix(data=c(1,1,10,10,100,100), nrow=3, ncol=2, byrow=TRUE))
#'
#' @seealso
#' \code{\link{Inf.Array}} for informative array testing
#'
#' \code{\link{NI.A2M}} for non-informative array testing with master pooling
#'
#' \code{\link{OTC}} for finding the optimal testing configuration for a number
#' of standard group testing algorithms
#'
#' \url{http://chrisbilder.com/grouptesting/HBTM}
#'
#' @references
#' \emph{Graff, L.E. & Roeloffs, R. (1972)}. Group testing in the presence of test
#' error; an extension of the Dorfman procedure. \emph{Technometrics, 14, 113-122}.
#'
#' \emph{Kim, H.Y.; Hudgens, M.G.; Dreyfuss, J.M.; Westreich, D.J. & Pilcher, C.D. (2007)}.
#' Comparison of group testing algorithms for case identification in the presence of test error.
#' \emph{Biometrics, 63, 1152-1163}.
#'
#' \emph{Malinovsky, Y.; Albert, P.S. & Roy, A. (2016)}. Reader reaction: A note on the
#' evaluation of group testing algorithms in the presence of misclassification.
#' \emph{Biometrics, 72, 299-302}.
#'
#' @family Optimal Testing Configuration functions

#    Brianna Hitt - 05-01-17

NI.Array <- function(p, Se, Sp, group.sz, obj.fn, weights){

  start.time<-base::proc.time()

  set.of.I <- group.sz

  save.it <- matrix(data=NA, nrow=length(set.of.I), ncol=18)
  count <- 1

  for(I in set.of.I){
    N <- I^2

    # build a matrix of probabilities
    # this is the same for an overall probability p and for a vector p
    p.mat <- matrix(data = p[1], nrow = I, ncol = I)

    # call Array.Measures to calculate descriptive measures for the given array size
    save.info <- Array.Measures(p = p.mat,se = Se, sp = Sp)

    # extract accuracy measures for each individual
    ET <- save.info$T
    PSe.mat <- save.info$PSe
    PSp.mat <- save.info$PSp
    MAR <- MAR.func(ET, PSe.mat, PSp.mat, p.mat)

    # calculate overall accuracy measures
    PSe <- sum(p.mat*PSe.mat)/sum(p.mat)
    PSp <- sum((1-p.mat)*(PSp.mat))/sum(1-p.mat)
    PPPV <- sum(p.mat*PSe.mat)/sum(p.mat*PSe.mat + (1-p.mat)*(1-PSp.mat))
    PNPV <- sum((1-p.mat)*PSp.mat)/sum((1-p.mat)*PSp.mat + p.mat*(1-PSe.mat))

    # for each row in the matrix of weights, calculate the GR function
    if(is.null(dim(weights))){
      GR1 <- NA
      GR2 <- NA
      GR3 <- NA
      GR4 <- NA
      GR5 <- NA
      GR6 <- NA
    } else{
      GR1 <- GR.func(ET, p.mat, PSe.mat, PSp.mat, D1=weights[1,1], D2=weights[1,2])
      if(dim(weights)[1]>=2){
        GR2 <- GR.func(ET, p.mat, PSe.mat, PSp.mat, D1=weights[2,1], D2=weights[2,2])
      } else{GR2 <- NA}
      if(dim(weights)[1]>=3){
        GR3 <- GR.func(ET, p.mat, PSe.mat, PSp.mat, D1=weights[3,1], D2=weights[3,2])
      } else{GR3 <- NA}
      if(dim(weights)[1]>=4){
        GR4 <- GR.func(ET, p.mat, PSe.mat, PSp.mat, D1=weights[4,1], D2=weights[4,2])
      } else{GR4 <- NA}
      if(dim(weights)[1]>=5){
        GR5 <- GR.func(ET, p.mat, PSe.mat, PSp.mat, D1=weights[5,1], D2=weights[5,2])
      } else{GR5 <- NA}
      if(dim(weights)[1]>=6){
        GR6 <- GR.func(ET, p.mat, PSe.mat, PSp.mat, D1=weights[6,1], D2=weights[6,2])
      } else{GR6 <- NA}
    }

    save.it[count,] <- c(p[1], Se, Sp, I, N, ET, ET/N, MAR, GR1/N, GR2/N, GR3/N, GR4/N, GR5/N, GR6/N, PSe, PSp, PPPV, PNPV)
    cat("Row/Column Size = ", I, ", Array Size = ", N, "\n", sep="")
    count <- count + 1

  }

  # find the optimal testing configuration, over all array sizes considered
  result.ET <- save.it[save.it[,7]==min(save.it[,7]),c(1:6,7,15:ncol(save.it))]
  result.MAR <- save.it[save.it[,8]==min(save.it[,8]),c(1:6,8,15:ncol(save.it))]
  result.GR1 <- save.it[save.it[,9]==min(save.it[,9]),c(1:6,9,15:ncol(save.it))]
  result.GR2 <- save.it[save.it[,10]==min(save.it[,10]),c(1:6,10,15:ncol(save.it))]
  result.GR3 <- save.it[save.it[,11]==min(save.it[,11]),c(1:6,11,15:ncol(save.it))]
  result.GR4 <- save.it[save.it[,12]==min(save.it[,12]),c(1:6,12,15:ncol(save.it))]
  result.GR5 <- save.it[save.it[,13]==min(save.it[,13]),c(1:6,13,15:ncol(save.it))]
  result.GR6 <- save.it[save.it[,14]==min(save.it[,14]),c(1:6,14,15:ncol(save.it))]

  p.mat.ET <- matrix(data = result.ET[1], nrow = result.ET[4], ncol = result.ET[4])
  p.mat.MAR <- matrix(data = result.MAR[1], nrow = result.MAR[4], ncol = result.MAR[4])
  if(is.null(dim(weights))){
    p.mat.GR1 <- NA
    p.mat.GR2 <- NA
    p.mat.GR3 <- NA
    p.mat.GR4 <- NA
    p.mat.GR5 <- NA
    p.mat.GR6 <- NA
  } else{
    p.mat.GR1 <- matrix(data = result.GR1[1], nrow = result.GR1[4], ncol = result.GR1[4])
    if(dim(weights)[1]>=2){
      p.mat.GR2 <- matrix(data = result.GR2[1], nrow = result.GR2[4], ncol = result.GR2[4])
    } else{p.mat.GR2 <- NA}
    if(dim(weights)[1]>=3){
      p.mat.GR3 <- matrix(data = result.GR3[1], nrow = result.GR3[4], ncol = result.GR3[4])
    } else{p.mat.GR3 <- NA}
    if(dim(weights)[1]>=4){
      p.mat.GR4 <- matrix(data = result.GR4[1], nrow = result.GR4[4], ncol = result.GR4[4])
    } else{p.mat.GR4 <- NA}
    if(dim(weights)[1]>=5){
      p.mat.GR5 <- matrix(data = result.GR5[1], nrow = result.GR5[4], ncol = result.GR5[4])
    } else{p.mat.GR5 <- NA}
    if(dim(weights)[1]>=6){
      p.mat.GR6 <- matrix(data = result.GR6[1], nrow = result.GR6[4], ncol = result.GR6[4])
    } else{p.mat.GR6 <- NA}
  }

  # create a list of results for each objective function
  opt.ET <- list("OTC" = list("Array.dim" = result.ET[4], "Array.sz" = result.ET[5]), "p.mat" = p.mat.ET, "ET" = result.ET[6], "value" = result.ET[7], "PSe" = result.ET[8], "PSp" = result.ET[9], "PPPV" = result.ET[10], "PNPV" = result.ET[11])
  opt.MAR <- list("OTC" = list("Array.dim" = result.MAR[4], "Array.sz" = result.MAR[5]), "p.mat" = p.mat.MAR, "ET" = result.MAR[6], "value" = result.MAR[7], "PSe" = result.MAR[8], "PSp" = result.MAR[9], "PPPV" = result.MAR[10], "PNPV" = result.MAR[11])
  opt.GR1 <- list("OTC" = list("Array.dim" = result.GR1[4], "Array.sz" = result.GR1[5]), "p.mat" = p.mat.GR1, "ET" = result.GR1[6], "value" = result.GR1[7], "PSe" = result.GR1[8], "PSp" = result.GR1[9], "PPPV" = result.GR1[10], "PNPV" = result.GR1[11])
  opt.GR2 <- list("OTC" = list("Array.dim" = result.GR2[4], "Array.sz" = result.GR2[5]), "p.mat" = p.mat.GR2, "ET" = result.GR2[6], "value" = result.GR2[7], "PSe" = result.GR2[8], "PSp" = result.GR2[9], "PPPV" = result.GR2[10], "PNPV" = result.GR2[11])
  opt.GR3 <- list("OTC" = list("Array.dim" = result.GR3[4], "Array.sz" = result.GR3[5]), "p.mat" = p.mat.GR3, "ET" = result.GR3[6], "value" = result.GR3[7], "PSe" = result.GR3[8], "PSp" = result.GR3[9], "PPPV" = result.GR3[10], "PNPV" = result.GR3[11])
  opt.GR4 <- list("OTC" = list("Array.dim" = result.GR4[4], "Array.sz" = result.GR4[5]), "p.mat" = p.mat.GR4, "ET" = result.GR4[6], "value" = result.GR4[7], "PSe" = result.GR4[8], "PSp" = result.GR4[9], "PPPV" = result.GR4[10], "PNPV" = result.GR4[11])
  opt.GR5 <- list("OTC" = list("Array.dim" = result.GR5[4], "Array.sz" = result.GR5[5]), "p.mat" = p.mat.GR5, "ET" = result.GR5[6], "value" = result.GR5[7], "PSe" = result.GR5[8], "PSp" = result.GR5[9], "PPPV" = result.GR5[10], "PNPV" = result.GR5[11])
  opt.GR6 <- list("OTC" = list("Array.dim" = result.GR6[4], "Array.sz" = result.GR6[5]), "p.mat" = p.mat.GR6, "ET" = result.GR6[6], "value" = result.GR6[7], "PSe" = result.GR6[8], "PSp" = result.GR6[9], "PPPV" = result.GR6[10], "PNPV" = result.GR6[11])

  # create a list of results, including all objective functions
  opt.all <- list("opt.ET" = opt.ET, "opt.MAR" = opt.MAR, "opt.GR1" = opt.GR1, "opt.GR2" = opt.GR2,
                  "opt.GR3" = opt.GR3, "opt.GR4" = opt.GR4, "opt.GR5" = opt.GR5, "opt.GR6" = opt.GR6)
  # remove any objective functions not requested by the user
  opt.req <- Filter(function(x) !is.na(x$ET), opt.all)
  time.it(start.time)
  c("prob" = list(p), "Se" = Se, "Sp" = Sp, opt.req)

}

###################################################################
