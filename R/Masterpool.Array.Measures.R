
# Start  MasterPool.Array.Measures() functions
###################################################################
#    Brianna Hitt - 05-01-17
#    Purpose: calculates descriptive measures for array testing with master pooling
#      calls: f, beta0, beta1, gamma0, gamma1, beta0.y, beta1.y, sum.beta.gamma,
#             mu.y, and nu.y - functions given in the equations for array testing
#             with master pooling in "Comparison of Group Testing Algorithms for Case 
#             Identification in the Presence of Test Error" by Kim et al. (2007)
#      inputs: results = an object containing results from Array.Measures()
#              n = size of a row/column in the square array
#              pmat = matrix of individual probabilities
#              Se = sensitivity of the diagnostic test
#              Sp = specificity of the diagnostic test
#      outputs: list of the expected number of tests (ET), and measures of testing accuracy,
#               including PSe, PSp, PPPV, and PNPV
#

f <- function(n, p, Se, Sp){
  (1 - Sp)*(1 - p)^n + Se*(1 - (1 - p)^n)
}

beta0 <- function(n, c, p){
  choose(n,c)*((1 - p)^(n^2 - n*c + c))*(1 - (1 - p)^(n-1))^c
}

beta1 <- function(n, c, p){
  choose(n,c)*((1 - p)^(n^2 - n*c))*(1 - (1 - p)^n)^c - beta0(n, c, p)
}

gamma0 <- function(n, c, Se, Sp){
  (1 - Sp)*((1 - Se)^c)*(Sp^(n-c))
}

gamma1 <- function(n, c, Se, Sp){
  Se*(Sp^(n-c))*(1 - Se)^c
}

beta0.y <- function(n, c, p){
  beta0(n, c, p)/(1 - p)
}

beta1.y <- function(n, c, p){
  choose((n-1),c)*((1 - (1 - p)^n)^c)*((1 - p)^(n^2 - n*c - 1)) + choose((n-1), (c-1))*((1 - (1 - p)^n)^(c-1))*((1 - p)^(n^2 - n*c))*(1 - (1 - p)^(n-1)) - beta0.y(n, c, p)
}

sum.beta.gamma <- function(n,p,Se,Sp){
  sum <- 0
  for(c in 1:n){
    sum <- sum + beta0.y(n,c,p)*gamma0(n,c,Se,Sp) + beta1.y(n,c,p)*gamma1(n,c,Se,Sp)
  }
  return(sum)
}

mu.y <- function(n, p, Se, Sp){
  f(n=(n^2 - 2*n + 1),p,Se,Sp)*((1 - Sp)*(1 - p)^(n-1))^2 + (Se^2)*(1 - (1 - p)^(n-1))*((1 - Sp)*(1 - p)^(n-1) + f(n=(n - 1),p,Se,Sp))
}

nu.y <- function(n, p, Se, Sp){
  (1 - Sp)*beta0.y(n,c=0,p)*gamma0(n,c=0,Se,Sp) + Se*sum.beta.gamma(n,p,Se,Sp)
}

MasterPool.Array.Measures <- function(results, n, pmat, Se, Sp){
  
  # extract the measures from the results for Array.Measures()
  ET.A2 <- results$T
  PSe.A2 <- results$PSe
  PSp.A2 <- results$PSp
  PPV.A2 <- results$PPV
  NPV.A2 <- results$NPV
  
  p <- pmat[1,1]
  
  # calculate the measures for array testing with master pooling
  ET.A2M <- 1/(n^2) + (1 - Se - Sp)*((1 - p)^(n^2))*(2/n + (1 - Sp)^2 + 2*(1 - Sp)*Sp^n) + Se*ET.A2
  PSe.A2M <- Se^4 + 2*(Se^3)*(1 - Se)*(1 - f(n,pmat,Se,Sp))^(n-1)
  PSp.A2M <- 1 - ((1 - Sp)*mu.y(n, pmat, Se, Sp) + 2*(1 - Sp)*nu.y(n, pmat, Se, Sp))
  PPV.A2M <- (pmat*PSe.A2M)/((1-pmat)*(1-PSp.A2M) + pmat*PSe.A2M)
  NPV.A2M <- ((1-pmat)*PSp.A2M)/(pmat*(1-PSe.A2M) + (1-pmat)*PSp.A2M)
  
  list("ET" = ET.A2M, "PSe" = PSe.A2M, "PSp" = PSp.A2M, "PPV" = PPV.A2M, "NPV" = NPV.A2M)
}

###################################################################
