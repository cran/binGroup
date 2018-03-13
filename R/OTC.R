#' @title Find the optimal testing configuration for standard group testing algorithms
#'
#' @author Brianna Hitt
#'
#' @description Find the optimal testing configuration (OTC) for standard group testing
#' algorithms, including hierarchical and array-based algorithms, and calculate the
#' associated operating characteristics.
#'
#' @details This function finds the optimal testing configuration and computes the
#' associated operating characteristics for standard group testing algorithms,
#' as described in Hitt et al. (2018) at \url{http://www.chrisbilder.com/grouptesting/HBTM/}.
#'
#' Available algorithms include two- and three-stage hierarchical testing, and
#' array testing with and without master pooling. Both non-informative and informative
#' group testing settings are allowed for each algorithm, with one exception. Only
#' non-informative array testing with master pooling is available, because no informative
#' group testing algorithms have been proposed for array testing with master pooling.
#' Operating characteristics calculated include the expected testing expenditure and
#' accuracy measures, including pooling sensitivity, pooling specificity, pooling
#' positive predictive value, and pooling negative predictive value for each individual.
#'
#' @param algorithm character string defining the group testing algorithm to be used.
#' Non-informative testing options include two-stage ("\kbd{D2}") and three-stage
#' ("\kbd{D3}") hierarchical, and square array testing with and without master pooling
#' ("\kbd{A2M}" and "\kbd{A2}", respectively). Informative testing options include
#' two-stage ("\kbd{ID2}") and three-stage ("\kbd{ID3}") hierarchical, and square
#' array testing without master pooling ("\kbd{IA2}").
#' @param p overall probability of disease that will be used to generate either a
#' homogeneous or heterogeneous vector/matrix of individual probabilities,
#' depending on the algorithm specified
#' @param probabilities a vector of individual probabilities, which is homogeneous for
#' non-informative testing algorithms and heterogeneous for informative testing algorithms
#' @param Se the sensitivity of the diagnostic test
#' @param Sp the specificity of the diagnostic test
#' @param group.sz a single group size for which to calculate operating characteristics
#' (for two-stage hierarchical and array testing) or to find the optimal testing
#' configuration over all possible configurations (for three-stage hierarchical testing),
#' or a range of group sizes over which to find the optimal testing configuration
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
#' @param alpha a shape parameter for the beta distribution that specifies the degree of
#' heterogeneity for the generated probability vector (for informative testing only)
#'
#' @note Either \kbd{p} or \kbd{probabilities} should be specified, but not both.
#' @return A list containing:
#' \item{prob}{the probability of infection, as specified by the user}
#' \item{alpha}{level of heterogeneity for the generated probability vector
#' (for informative testing only)}
#' \item{Se}{the sensitivity of the diagnostic test}
#' \item{Sp}{the specificity of the diagnostic test}
#' \item{opt.ET, opt.MAR, opt.GR}{a list for each objective function specified
#' by the user, containing:
#' \describe{
#' \item{OTC}{a list specifying the optimal testing configuration, which may include:
#' \describe{
#' \item{Stage1}{pool sizes for the first stage of testing, if applicable}
#' \item{Stage2}{pool sizes for the second stage of testing, if applicable}
#' \item{Block.sz}{the block size/initial group size for informative Dorfman testing,
#' which is not tested}
#' \item{pool.szs}{pool sizes for the first stage of testing for informative Dorfman
#' testing}
#' \item{Array.dim}{the row/column size for array testing}
#' \item{Array.sz}{the array size for array testing}}}
#' \item{p.vec}{the sorted vector of individual probabilities, if applicable}
#' \item{p.mat}{the sorted matrix of individual probabilities in gradient arrangement,
#' if applicable}
#' \item{ET}{the expected testing expenditure for the OTC}
#' \item{value}{the value of the objective function per individual}
#' \item{PSe}{the overall pooling sensitivity for the algorithm}
#' \item{PSp}{the overall pooling specificity for the algorithm}
#' \item{PPPV}{the overall pooling positive predictive value for the algorithm}
#' \item{PNPV}{the overall pooling negative predictive value for the algorithm}}}
#'
#' @examples
#' # Find the optimal testing configuration for non-informative
#' #   two-stage hierarchical (Dorfman) testing
#' # This example takes less than one second to run
#' OTC(algorithm="D2", p=0.05, Se=0.99, Sp=0.99, group.sz=2:100,
#' obj.fn=c("ET", "MAR"))
#'
#' # Find the optimal testing configuration for informative
#' #   two-stage hierarchical (Dorfman) testing, implemented
#' #   via the pool-specific optimal Dorfman (PSOD) method
#' #   described in McMahan et al. (2012a)
#' # This example takes approximately 16 minutes to run
#' \dontrun{set.seed(52613)
#' OTC(algorithm="ID2", p=0.01, Se=0.95, Sp=0.95, group.sz=3:50,
#' obj.fn=c("ET", "MAR", "GR"),
#' weights=matrix(data=c(1, 1, 10, 10, 0.5, 0.5), nrow=3, ncol=2, byrow=TRUE),
#' alpha=0.5)}
#'
#' # Find the optimal testing configuration over all possible
#' #   testing configurations for a specified group size for
#' #   non-informative three-stage hierarchical testing
#' # This example takes approximately 1 second to run
#' OTC(algorithm="D3", p=0.001, Se=0.95, Sp=0.95, group.sz=24,
#' obj.fn=c("ET", "MAR", "GR"),
#' weights=matrix(data=c(1, 1), nrow=1, ncol=2, byrow=TRUE))
#'
#' # Find the optimal testing configuration for non-informative
#' #   three-stage hierarchical testing
#' # This example takes approximately 20 seconds to run
#' \dontrun{OTC(algorithm="D3", p=0.06, Se=0.90, Sp=0.90, group.sz=3:30,
#' obj.fn=c("ET", "MAR", "GR"),
#' weights=matrix(data=c(1, 1, 10, 10, 100, 100), nrow=3, ncol=2, byrow=TRUE))}
#'
#' # Find the optimal testing configuration over all possible
#' #   configurations for a specified group size, given a heterogeneous
#' #   vector of probabilities
#' # This example takes less than 1 second to run
#' OTC(algorithm="ID3",
#' probabilities=c(0.012, 0.014, 0.011, 0.012, 0.010, 0.015),
#' Se=0.99, Sp=0.99, group.sz=6, obj.fn=c("ET","MAR","GR1"),
#' weights=matrix(data=c(1, 1), nrow=1, ncol=2, byrow=TRUE), alpha=0.5)
#'
#' # Calculate the operating characteristics for a specified array size
#' #   for non-informative array testing without master pooling
#' # This example takes less than 1 second to run
#' OTC(algorithm="A2", p=0.005, Se=0.95, Sp=0.95, group.sz=15, obj.fn=c("ET", "MAR"))
#'
#' # Find the optimal testing configuration for informative
#' #   array testing without master pooling
#' # This example takes approximately 30 seconds to run
#' \dontrun{set.seed(1002)
#' OTC(algorithm="IA2", p=0.03, Se=0.95, Sp=0.95, group.sz=3:20,
#' obj.fn=c("ET", "MAR", "GR"),
#' weights=matrix(data=c(1, 1, 10, 10, 100, 100), nrow=3, ncol=2, byrow=TRUE),
#' alpha=2)}
#'
#' # Find the optimal testing configuration for non-informative
#' #   array testing with master pooling
#' # This example takes approximately 15 seconds to run
#' \dontrun{OTC(algorithm="A2M", p=0.02, Se=0.90, Sp=0.90, group.sz=3:20, obj.fn=c("ET", "MAR", "GR"),
#' weights=matrix(data=c(1, 1, 10, 10, 0.5, 0.5, 2, 2, 100, 100, 10, 100), nrow=6, ncol=2, byrow=TRUE))}
#'
#' @seealso
#' \code{\link{NI.Dorf}} for non-informative two-stage (Dorfman) testing, \code{\link{Inf.Dorf}} for
#' informative two-stage (Dorfman) testing, \code{\link{NI.D3}} for non-informative three-stage
#' hierarchical testing, \code{\link{Inf.D3}} for informative three-stage hierarchical testing,
#' \code{\link{NI.Array}} for non-informative array testing, \code{\link{Inf.Array}} for informative
#' array testing, and \code{\link{NI.A2M}} for non-informative array testing with master pooling
#'
#' \url{http://chrisbilder.com/grouptesting/HBTM/}
#'
#' @references
#' \emph{Graff, L.E. & Roeloffs, R. (1972)}. Group testing in the presence of test error; an extension of
#' the Dorfman procedure. \emph{Technometrics, 14, 113-122}.
#'
#' \emph{Malinovsky, Y.; Albert, P.S. & Roy, A. (2016)}. Reader reaction: A note on the evaluation of
#' group testing algorithms in the presence of misclassification. \emph{Biometrics, 72, 299-302}.
#'
#' @family Optimal Testing Configuration functions

OTC <- function(algorithm, p=NULL, probabilities=NULL, Se=0.99, Sp=0.99, group.sz, obj.fn=c("ET","MAR"), weights=NULL, alpha=NULL){

  ## make sure that all necessary information is included in the correct format
  if(is.null(p) & is.null(probabilities)){
    stop("Please specify an overall probability of disease using the 'p' argument, \n or specify a vector of individual probabilities using the 'probabilities' argument.")
  } else if(!is.null(p) & !is.null(probabilities)){
    stop("You have specified both an overall probability of disease AND a \n vector of individual probabilities. Please specify only one option.")
  } else{
    if(!is.null(p)){
      if(length(p)==1){
        cat("You have specified an overall probability of disease. \n A probability vector will be generated based on the algorithm specified.\n")
      } else{
        stop("You have specified a probability vector instead of an overall probability of disease.\n Please specify an overall probability of disease, and the probability vector will be \n generated based on the algorithm specified for each group size included in the range.\n")
      }
      if((algorithm %in% c("ID2", "ID3", "IA2")) & is.null(alpha)){
        stop("Please specify the level of heterogeneity for generating the vector of individual probabilities using the 'alpha' argument.")
      }
    }
    if(!is.null(probabilities)){
      if(length(group.sz)==1){
        cat("You have specified a vector containing individual probabilities of disease.\n")
        if((algorithm %in% c("D2","D3","ID2","ID3")) & length(probabilities)!=group.sz){
          stop("The vector of individual probabilities is not the correct length. Please make sure\n that the length of the probability vector is the same as the specified group size.\n")
        } else if((algorithm %in% c("A2","A2M","IA2")) & length(probabilities)!=group.sz^2){
          stop("The vector of individual probabilities is not the correct length. Please make sure that the\n length of the probability vector is the same as the overall array size (the specified group size squared).\n")
        }
        if((algorithm %in% c("D2","D3","A2","A2M")) & all.equal(probabilities, rep(probabilities[1],length(probabilities)))!=TRUE){
          stop("You have specified a heterogeneous probability vector for a non-informative\n algorithm. Please specify a homogeneous probability vector using the 'probabilities'\n argument or specify an overall probability of disease using the 'p' argument.\n")
        }
      } else if(length(group.sz)>1){
        stop("You have specified a probability vector along with a range \n of group sizes. Please specify a single group size.\n")
      }
      if(!is.null(alpha)){
        cat("You have specified a vector of individual probabilities - alpha will be ignored.")
      }
    }
  }

  if(length(group.sz)==1){
    if(algorithm %in% c("D3","ID2")){
      cat("A single group size was provided. The optimal testing configuration will be found \n over all possible testing configurations for the specified group size.\n")
    } else{
      cat("A single group size was provided. No optimization will be performed.\n")
    }
  }

  if(is.null(obj.fn)){
    stop("Please specify one or more objective functions for which to find the optimal testing configuration.\n")
  }

  if("GR" %in% obj.fn){
    if(is.null(weights)){
      stop("No weights have been specified. The GR function will not be calculated.\n")
    } else if(dim(weights)[2]!=2){
      stop("Please check the dimension of the weights matrix. \n Each row should specify a set of weights, D1 and D2.\n")
    }
  }

  # call function for non-informative two-stage hierarchical (Dorfman) testing
  if(algorithm == "D2"){
    if(min(group.sz)<2){
      stop("Please specify a minimum group size of at least 2.\n")
    }
    if(!is.null(alpha)){
      cat("Non-informative Testing - alpha will be ignored. \n")
    }
    cat("Algorithm: Non-informative two-stage hierarchical (Dorfman) testing \n")
    if(!is.null(p)){
      results <- NI.Dorf(p, Se, Sp, group.sz, obj.fn, weights)
    } else if(!is.null(probabilities)){
      results <- NI.Dorf(probabilities, Se, Sp, group.sz, obj.fn, weights)
    }
  }

  # call function for non-informative three-stage hierarchical testing
  if(algorithm == "D3"){
    if(min(group.sz)<3){
      stop("Please specify a minimum group size of at least 3.\n")
    }
    if(max(group.sz)>50){
      message("NOTE: You have specified a maximum group size larger than 50.\n This function may take a VERY long time to run.\n Press 'ESC' if you wish to cancel the submitted statements.\n")
    }
    if(!is.null(alpha)){
      cat("Non-informative Testing - alpha will be ignored. \n")
    }
    cat("Algorithm: Non-informative three-stage hierarchical testing \n")
    if(!is.null(p)){
      results <- NI.D3(p, Se, Sp, group.sz, obj.fn, weights)
    } else if(!is.null(probabilities)){
      results <- NI.D3(probabilities, Se, Sp, group.sz, obj.fn, weights)
    }
  }

  # call function for non-informative square array testing without master pooling
  if(algorithm == "A2"){
    if(min(group.sz)<3){
      stop("Please specify a minimum group size of at least 3.\n")
    }
    if(max(group.sz)>50){
      message("NOTE: You have specified a maximum group size larger than 50.\n This function may take a VERY long time to run.\n Press 'ESC' if you wish to cancel the submitted statements.\n")
    }
    if(!is.null(alpha)){
      cat("Non-informative Testing - alpha will be ignored. \n")
    }
    cat("Algorithm: Non-informative square array testing without master pooling \n")
    if(!is.null(p)){
      results <- NI.Array(p, Se, Sp, group.sz, obj.fn, weights)
    } else if(!is.null(probabilities)){
      results <- NI.Array(probabilities, Se, Sp, group.sz, obj.fn, weights)
    }
  }

  # call function for non-informative square array testing with master pooling
  if(algorithm == "A2M"){
    if(min(group.sz)<3){
      stop("Please specify a minimum group size of at least 3.\n")
    }
    if(max(group.sz)>50){
      message("NOTE: You have specified a maximum group size larger than 50.\n This function may take a VERY long time to run.\n Press 'ESC' if you wish to cancel the submitted statements.\n")
    }
    if(!is.null(alpha)){
      cat("Non-informative Testing - alpha will be ignored. \n")
    }
    cat("Algorithm: Non-informative square array testing with master pooling \n")
    if(!is.null(p)){
      results <- NI.A2M(p, Se, Sp, group.sz, obj.fn, weights)
    } else if(!is.null(probabilities)){
      results <- NI.A2M(probabilities, Se, Sp, group.sz, obj.fn, weights)
    }
  }

  # call function for informative two-stage hierarchical (Dorfman) testing
  if(algorithm == "ID2"){
    if(min(group.sz)<2){
      stop("Please specify a minimum group size of at least 2.\n")
    }
    if(max(group.sz)>=50){
      message("NOTE: You have specified a maximum group size larger than 50.\n This function may take a VERY long time to run.\n Press 'ESC' if you wish to cancel the submitted statements.\n")
    }
    cat("Algorithm: Informative Dorfman testing \n")
    if(!is.null(p)){
      results <- Inf.Dorf(p, Se, Sp, group.sz, obj.fn, weights, alpha)
    } else if(!is.null(probabilities)){
      results <- Inf.Dorf(probabilities, Se, Sp, group.sz, obj.fn, weights, alpha=NA)
    }
  }

  # call function for informative three-stage hierarchical testing
  if(algorithm == "ID3"){
    if(min(group.sz)<3){
      stop("Please specify a minimum group size of at least 3.\n")
    }
    if(max(group.sz)>=50){
      message("NOTE: You have specified a maximum group size larger than 50.\n This function may take a VERY long time to run.\n Press 'ESC' if you wish to cancel the submitted statements.\n")
    }
    cat("Algorithm: Informative three-stage hierarchical testing \n")
    if(!is.null(p)){
      results <- Inf.D3(p, Se, Sp, group.sz, obj.fn, weights, alpha)
    } else if(!is.null(probabilities)){
      results <- Inf.D3(probabilities, Se, Sp, group.sz, obj.fn, weights, alpha=NA)
    }
  }

  # call function for informative square array testing without master pooling
  if(algorithm == "IA2"){
    if(min(group.sz)<3){
      stop("Please specify a minimum group size of at least 3.\n")
    }
    if(max(group.sz)>=50){
      message("NOTE: You have specified a maximum group size larger than 50.\n This function may take a VERY long time to run.\n Press 'ESC' if you wish to cancel the submitted statements.\n")
    }
    cat("Algorithm: Informative square array testing without master pooling \n")
    if(!is.null(p)){
      results <- Inf.Array(p, Se, Sp, group.sz, obj.fn, weights, alpha)
    } else if(!is.null(probabilities)){
      results <- Inf.Array(probabilities, Se, Sp, group.sz, obj.fn, weights, alpha=NA)
    }
  }

  results
}



# Supporting functions for OTC and the associated calls
####################################################################

p.vec.func <- function(p, alpha, grp.sz){
  if(is.na(p)){
    NA
  } else{
    p.try <- try(beta.dist(p=p, alpha=alpha, grp.sz=grp.sz), silent=TRUE)
    if(class(p.try)=="try-error"){
      beta.dist(p=p, alpha=alpha, grp.sz=grp.sz, simul=TRUE)
    } else{
      beta.dist(p=p, alpha=alpha, grp.sz=grp.sz)
    }
  }
}

###################################################################




# Start MAR.func() function
###################################################################
#    Brianna Hitt - 4-17-17
#    Purpose: calculates MAR objective function, from Malinovsky, Albert & Roy (2015)
#      inputs: ET - expected number of tests
#              p.vec - vector of individual probabilities
#              PSe.vec - vector of individual pooling sensitivities
#              PSp.vec - vector of individual pooling specificities
#      note: The MAR objective function divides ET, the expected number of tests, by EC,
#            the expected number of correct classifications, and should be minimized.
#            Note: Malinovsky, Albert, & Roy (2015) maximized the reciprocal, E(C)/E(T).

MAR.func <- function(ET, p.vec, PSe.vec, PSp.vec){
  EC <- sum(PSe.vec*p.vec + PSp.vec*(1 - p.vec))
  ET/EC
}
###################################################################




# Start GR.func() function
###################################################################
#    Brianna Hitt - 4-17-17
#    Purpose: calculates GR objective function, from Graff & Roeloffs (1972)
#               M = E(T) + D_1*(# of misclassified negatives) + D_2*(# of misclassified positives)
#      inputs: ET - expected number of tests
#              p.vec - vector of individual probabilities
#              PSe.vec - vector of individual pooling sensitivities
#              PSp.vec - vector of individual pooling specificities
#              D1, D2 - weights/costs for misclassification
#      note: this function specifies equal weights of 1 by default

GR.func <- function(ET, p.vec, PSe.vec, PSp.vec, D1=1, D2=1){
  ET + D1*sum((1-PSp.vec)*(1-p.vec)) + D2*sum((1-PSe.vec)*p.vec)
}
###################################################################




# Start time.it() function
###################################################################
#    Brianna Hitt - 5-13-17
#    Purpose: calculates the time elapsed
#      inputs: x = object containing the start time

time.it <- function(x) {
  end.time<-base::proc.time()
  save.time<-end.time-x
  cat("\n Number of minutes running:", save.time[3]/60, "\n \n")
  save.time[3]/60
}
###################################################################

