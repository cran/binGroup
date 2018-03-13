
# Start  NI.Dorf() function
###################################################################

#' @title Non-informative two-stage hierarchical testing
#'
#' @author Brianna Hitt
#'
#' @description Find the optimal testing configuration for non-informative
#' two-stage hierarchical (Dorfman) testing and calculate the associated
#' operating characteristics.
#'
#' @details This function finds the optimal testing configuration and computes
#' the associated operating characteristics for non-informative two-stage
#' hierarchical (Dorfman) testing. See Hitt et al. (2018) at
#' \url{http://www.chrisbilder.com/grouptesting/HBTM/}, Dorfman (1943), or
#' Kim et al. (2007) for additional details.
#'
#' @param p the probability of disease, which can be specified as an overall
#' probability of disease or a homogeneous vector of individual probabilities
#' @param Se the sensitivity of the diagnostic test
#' @param Sp the specificity of the diagnostic test
#' @param group.sz a single group size for which to calculate the operating
#' characteristics, or a range of group sizes over which to find the optimal
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
#'
#' @return A list containing:
#' \item{prob}{the probability of disease, as specified by the user}
#' \item{Se}{the sensitivity of the diagnostic test}
#' \item{Sp}{the specificity of the diagnostic test}
#' \item{opt.ET, opt.MAR, opt.GR}{a list for each objective function specified by the
#' user, containing:
#' \describe{
#' \item{OTC}{a list specifying the optimal testing configuration, which includes:
#' \describe{
#' \item{Stage1}{the pool size for the first stage of testing, i.e. the initial
#' group size}}}
#' \item{p.vec}{the vector of individual probabilities}
#' \item{ET}{the expected testing expenditure for the OTC}
#' \item{value}{the value of the objective function per individual}
#' \item{PSe}{the overall pooling sensitivity for the algorithm}
#' \item{PSp}{the overall pooling specificity for the algorithm}
#' \item{PPPV}{the overall pooling positive predictive value for the algorithm}
#' \item{PNPV}{the overall pooling negative predictive value for the algorithm}}}
#'
#' @examples
#' # Find the optimal testing configuration for non-informative two-stage
#' #   hierarchical (Dorfman) testing over a range of group sizes
#' NI.Dorf(p=0.01, Se=0.95, Sp=0.95, group.sz=2:100, obj.fn=c("ET", "MAR"),
#' weights=NULL)
#'
#' # Calculate the operating characteristics for a specified initial
#' #   group size for non-informative two-stage hierarchical (Dorfman) testing
#' NI.Dorf(p=rep(0.025, 50), Se=0.90, Sp=0.90, group.sz=50,
#' obj.fn=c("ET", "MAR", "GR"),
#' weights=matrix(data=c(1,1,10,10), nrow=2, ncol=2, byrow=TRUE))
#'
#' @seealso
#' \code{\link{Inf.Dorf}} for informative two-stage hierarchical (Dorfman) testing
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
#' \emph{Kim, H.Y.; Hudgens, M.G.; Dreyfuss, J.M.; Westreich, D.J. & Pilcher, C.D. (2007)}.
#' Comparison of group testing algorithms for case identification in the presence of test error.
#' \emph{Biometrics, 63, 1152-1163}.
#'
#' \emph{Malinovsky, Y.; Albert, P.S. & Roy, A. (2016)}. Reader reaction: A note on the
#' evaluation of group testing algorithms in the presence of misclassification.
#' \emph{Biometrics, 72, 299-302}.
#'
#' @family Optimal Testing Configuration functions

#    Brianna Hitt - 4-17-17
NI.Dorf <- function(p, Se, Sp, group.sz, obj.fn, weights){

  start.time<-base::proc.time()

  set.of.I <- group.sz
  save.it <- matrix(data = NA, nrow = length(set.of.I), ncol = 17)
  count <- 1

  for(I in set.of.I){
    # generate a probability vector for homogeneous population
    p.vec <- rep(x = p[1], times = I)

    # calculate descriptive measures for two-stage hierarchical testing
    save.info <- hierarchical.desc2(p = p.vec, se = Se, sp = Sp, I2 = NULL, order.p = FALSE)

    # extract ET, PSe, PSp and calculate the MAR function
    ET <- save.info$ET
    PSe.vec <- save.info$individual.testerror$pse.vec
    PSp.vec <- save.info$individual.testerror$psp.vec
    MAR <- MAR.func(ET, p.vec, PSe.vec, PSp.vec)

    # for non-informative Dorfman (two-stage hierarchical) testing, all individuals have the same testing accuracy measures
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

    save.it[count,] <- c(p[1], Se, Sp, I, ET, ET/I, MAR, GR1/I, GR2/I, GR3/I, GR4/I, GR5/I, GR6/I, PSe, PSp, PPPV, PNPV)
    cat("Initial Group Size =", I, "\n")
    count <- count + 1
  }

  # find the testing configuration with the minimum value, for each objective function
  result.ET <- save.it[save.it[,6]==min(save.it[,6]),c(1:5,6,14:ncol(save.it))]
  result.MAR <- save.it[save.it[,7]==min(save.it[,7]),c(1:5,7,14:ncol(save.it))]
  result.GR1 <- save.it[save.it[,8]==min(save.it[,8]),c(1:5,8,14:ncol(save.it))]
  result.GR2 <- save.it[save.it[,9]==min(save.it[,9]),c(1:5,9,14:ncol(save.it))]
  result.GR3 <- save.it[save.it[,10]==min(save.it[,10]),c(1:5,10,14:ncol(save.it))]
  result.GR4 <- save.it[save.it[,11]==min(save.it[,11]),c(1:5,11,14:ncol(save.it))]
  result.GR5 <- save.it[save.it[,12]==min(save.it[,12]),c(1:5,12,14:ncol(save.it))]
  result.GR6 <- save.it[save.it[,13]==min(save.it[,13]),c(1:5,13,14:ncol(save.it))]

  p.vec.ET <- rep(x = result.ET[1], times = result.ET[4])
  p.vec.MAR <- rep(x = result.MAR[1], times = result.MAR[4])
  if(is.null(dim(weights))){
    p.vec.GR1 <- NA
    p.vec.GR2 <- NA
    p.vec.GR3 <- NA
    p.vec.GR4 <- NA
    p.vec.GR5 <- NA
    p.vec.GR6 <- NA
  } else{
    p.vec.GR1 <- rep(x = result.GR1[1], times = result.GR1[4])
    if(dim(weights)[1]>=2){
      p.vec.GR2 <- rep(x = result.GR2[1], times = result.GR2[4])
    } else{p.vec.GR2 <- NA}
    if(dim(weights)[1]>=3){
      p.vec.GR3 <- rep(x = result.GR3[1], times = result.GR3[4])
    } else{p.vec.GR3 <- NA}
    if(dim(weights)[1]>=4){
      p.vec.GR4 <- rep(x = result.GR4[1], times = result.GR4[4])
    } else{p.vec.GR4 <- NA}
    if(dim(weights)[1]>=5){
      p.vec.GR5 <- rep(x = result.GR5[1], times = result.GR5[4])
    } else{p.vec.GR5 <- NA}
    if(dim(weights)[1]>=6){
      p.vec.GR6 <- rep(x = result.GR6[1], times = result.GR6[4])
    } else{p.vec.GR6 <- NA}
  }

  # create a list of results for each objective function
  opt.ET <- list("OTC" = list("Stage1" = result.ET[4]), "p.vec" = p.vec.ET, "ET" = result.ET[5], "value" = result.ET[6], "PSe" = result.ET[7], "PSp" = result.ET[8], "PPPV" = result.ET[9], "PNPV" = result.ET[10])
  opt.MAR <- list("OTC" = list("Stage1" = result.MAR[4]), "p.vec" = p.vec.MAR, "ET" = result.MAR[5], "value" = result.MAR[6], "PSe" = result.MAR[7], "PSp" = result.MAR[8], "PPPV" = result.MAR[9], "PNPV" = result.MAR[10])
  opt.GR1 <- list("OTC" = list("Stage1" = result.GR1[4]), "p.vec" = p.vec.GR1, "ET" = result.GR1[5], "value" = result.GR1[6], "PSe" = result.GR1[7], "PSp" = result.GR1[8], "PPPV" = result.GR1[9], "PNPV" = result.GR1[10])
  opt.GR2 <- list("OTC" = list("Stage1" = result.GR2[4]), "p.vec" = p.vec.GR2, "ET" = result.GR2[5], "value" = result.GR2[6], "PSe" = result.GR2[7], "PSp" = result.GR2[8], "PPPV" = result.GR2[9], "PNPV" = result.GR2[10])
  opt.GR3 <- list("OTC" = list("Stage1" = result.GR3[4]), "p.vec" = p.vec.GR3, "ET" = result.GR3[5], "value" = result.GR3[6], "PSe" = result.GR3[7], "PSp" = result.GR3[8], "PPPV" = result.GR3[9], "PNPV" = result.GR3[10])
  opt.GR4 <- list("OTC" = list("Stage1" = result.GR4[4]), "p.vec" = p.vec.GR4, "ET" = result.GR4[5], "value" = result.GR4[6], "PSe" = result.GR4[7], "PSp" = result.GR4[8], "PPPV" = result.GR4[9], "PNPV" = result.GR4[10])
  opt.GR5 <- list("OTC" = list("Stage1" = result.GR5[4]), "p.vec" = p.vec.GR5, "ET" = result.GR5[5], "value" = result.GR5[6], "PSe" = result.GR5[7], "PSp" = result.GR5[8], "PPPV" = result.GR5[9], "PNPV" = result.GR5[10])
  opt.GR6 <- list("OTC" = list("Stage1" = result.GR6[4]), "p.vec" = p.vec.GR6, "ET" = result.GR6[5], "value" = result.GR6[6], "PSe" = result.GR6[7], "PSp" = result.GR6[8], "PPPV" = result.GR6[9], "PNPV" = result.GR6[10])

  # create a list of results, including all objective functions
  opt.all <- list("opt.ET" = opt.ET, "opt.MAR" = opt.MAR, "opt.GR1" = opt.GR1, "opt.GR2" = opt.GR2,
                  "opt.GR3" = opt.GR3, "opt.GR4" = opt.GR4, "opt.GR5" = opt.GR5, "opt.GR6" = opt.GR6)
  # remove any objective functions not requested by the user
  opt.req <- Filter(function(x) !is.na(x$ET), opt.all)
  time.it(start.time)
  c("prob" = list(p), "Se" = Se, "Sp" = Sp, opt.req)

}

###################################################################
