\name{binGroup-package}
\Rdversion{1.1}
\alias{binGroup-package}
\alias{binGroup}
\docType{package}
\title{
Statistical Methods for Group Testing.
}
\description{
This package provides methods for estimation and hypothesis testing of proportions in group testing designs.
Currently, methods of estimating a proportion in a single population, as well as hypothesis tests and functions for experimental design for this situation are provided.
Further, regression methods are implemented for simple pooling and matrix pooling designs.
}
\details{
\tabular{ll}{
Package: \tab binGroup\cr
Type: \tab Package\cr
Version: \tab 1.0-3\cr
Date: \tab 2009-10-27\cr
License: \tab GPL\cr
LazyLoad: \tab no\cr
}

1) One-sample case

Methods for calculating confidence intervals for a single population proportion from designs with equal group sizes (as described by Tebbs and Bilder, 2004 and Schaarschmidt, 2007) are implemented in the function \code{\link{bgtCI}}.

If a confidence interval for a single proportion shall be calculated based on a design with a limited number of different group sizes, the exact method described by Hepworth (2004) is implemented in the function \code{\link{bgtvs}}.

For the problem of choosing an adequate experimental design in the one-sample case with only one group size, the functions \code{\link{estDesign}}, \code{\link{sDesign}}, \code{\link{nDesign}} implement different iterative approaches, as examplified by Swallow (1985) and Schaarschmidt (2007).

2) Regression models

Two approaches (by Vansteelandt et al., 2000 and Xie, 2001) to estimate parameters of group testing regression models can be applied by calling \code{\link{gtreg}}. Once fitted, corresponding methods to extract residuals, calculate predictions and summarize the parameter estimates (including hypotheses tests) are available in the S3 methods \code{\link{residuals.gt}}, \code{\link{predict.gt}} and \code{\link{summary.gt}}.

Group testing regression models in settings with matrix pooling (Xie, 2001) can be fit using \code{\link{gtreg.mp}}. 


}
\author{
Boan Zhang, Christopher Bilder, Frank Schaarschmidt

Maintainer: Frank Schaarschmidt <schaarschmidt@biostat.uni-hannover.de>
}

\references{

\emph{Hepworth, G. (1996)} Exact confidence intervals for proportions estimated by group testing. Biometrics 52, 1134-1146.

\emph{Schaarschmidt, F. (2007)} Experimental design for one-sided confidence intervals or hypothesis tests in binomial group testing. Communications in Biometry and Crop Science 2 (1), 32-40.
 http://agrobiol.sggw.waw.pl/cbcs/

\emph{Swallow, W.H. (1985)} Group testing for estimating infection rates and probabilities of disease transmission. Phytopathology Vol.75, N.8, 882-889.

\emph{Tebbs, J.M. & Bilder, C.R. (2004)} Confidence interval procedures for the probability of disease transmission in multiple-vector-transfer designs. Journal of Agricultural, Biological and Environmental Statistics, Vol.9, N.1, 75-90.

\emph{Vansteelandt, S., Goetghebeur, E., and Verstraeten, T. (2000)} Regression models for disease prevalence with diagnostic tests on pools of serum samples, \emph{Biometrics, 56, 1126-1133.}

\emph{Xie, M. (2001)} Regression analysis of group testing samples, \emph{Statistics in Medicine, 20, 1957-1969.}

}

\keyword{package}
\keyword{regression}
\keyword{models}
\keyword{htest}

\examples{

# 1) One-sample problem

# 1.1) Confidence intervals for designs with equal group size (pool size),
# where 
#  n denotes the number of groups (pools),
#  s denotes the common group size (number of individuals pooled per group),
#  Y denotes the number of groups tested positive.

# The following call reproduces the example given 
# by Tebbs and Bilder (2004) for the two-sided 95-percent 
# exact (Clopper-Pearson) interval:

bgtCI(n=24, Y=3, s=7, conf.level=0.95,
 alternative="two.sided", method="CP")

# 1.2) Confidence intervals for designs with unequal group size (pool size):
# Keeping notation as above but allowing for (a limited number of) different
# group size s, the examples given in Hepworth (1996), Table 5 can be 
# reproduced by calling:

 bgtvs(n=c(2,3), s=c(5,2), Y=c(0,0))
 bgtvs(n=c(2,3), s=c(5,2), Y=c(0,1))

# 1.3) For experimental design based on the bias of the point estimate,
# as proposed by Swallow (1985): The values in Table 1 (Swallow, 1985),
# p.885 can be reproduced by calling:

estDesign(n=10, smax=100, p.tr=0.001)
estDesign(n=10, smax=100, p.tr=0.01)


# 2.1) Fitting a regression model
# A HIV surveillance data (used by Vansteelandt et al. 2000)
# can be analysed for the dependence of HIV prevalence
# on covariates AGE and EDUC., with sensitivity and specificity
# assumed to be 0.9 each.

data(hivsurv)
fit1 <- gtreg(formula = groupres ~ AGE + EDUC., data = hivsurv,
  groupn = gnum, sens = 0.9, spec = 0.9, method = "Xie")
summary(fit1)


# 2.2) Fitting a regression model for matrix pooling data
# The function sim.mp is used to simulate a matrix pooling data set:

set.seed(9128)
sa1a<-sim.mp(beta.par=c(-7,0.1), rown=c(5,4), coln=c(6,5),
 sens=0.95, spec=0.95)
sa1<-sa1a$dframe

\dontrun{
fit1mp <- gtreg.mp(formula = cbind(col.resp, row.resp) ~ x, data = sa1, 
  col.groupn = coln, row.groupn = rown, arrayn = sqn, 
  sens = 0.95, spec = 0.95, linkf = "logit", n.gibbs = 1000, tol = 0.005)
summary(fit1mp)
}

}
