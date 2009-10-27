\name{bgtCI}

\alias{bgtCI}
\alias{bgtCP}
\alias{bgtBlaker}
\alias{bgtAC}
\alias{bgtWilson}
\alias{bgtWald}
\alias{bgtSOC}

\title{ Point estimator and Confidence intervals for a binomial proportion in binomial group testing  }

\description{ Calculates the point estimate, the exact Clopper-Pearson and Blaker CI, the Score test derived 
		Wilson and Agresti-Coull CI, the asymptotic second-order corrected interval fo Cai and the Wald CI for
		a single binomial proportion estimated from a binomial group testing trial. Assumes equal group sizes, 
		an assay method classifying a group as positive if at least one unit in the group is positive, individuals
		units randomly assigned to the groups.
 }
\usage{
bgtCI(n, s, Y, conf.level = 0.95,
 alternative = "two.sided", method = "CP")
}

\arguments{
  \item{n}{ integer, specifying the number of  groups (i.e. assays i.e. observations) }
  \item{s}{ integer, specifying the common size of groups i.e. the number of individual units in each group }
  \item{Y}{ integer, specifying the number of positive groups }
  \item{conf.level}{ nominal confidence level of the interval }
  \item{alternative}{ character string defining the alternative hypothesis, either 'two.sided', 'less' or 'greater'
  where 'less' gives the only an upper bound with confidence level=conf.level
  'greater' gives the only a lower bound with confidence level=conf.level
   and 'two.sided' gives a two-sided confidence interval with confidence level=conf.level}
   \item{method}{ character string defining the method for CI calculation, where: 
     \kbd{"CP"} is Clopper-Pearson, an exact tail interval showing symmetric coverage probability (inversion of two one-sided tests),
     \kbd{"Blaker"} is the Blaker interval, an exact interval, inversion of one two.sided test, therefore defined only two.sided, but shorter
             than the two-sided Clopper-Pearson CI. Both guarantee to contain the true parameter with at least conf.level*100 percent
             probability,
     \kbd{"AC"} is the Agresti-Coull (generalized Agresti-Coull) interval, asymptotic method,
     \kbd{"Score"} is Wilson Score, asymptotic method derived from inversion of the Score test,
     \kbd{"SOC"} is the second order corrected interval, asymptotic method for one-sided problems (for details see Cai, 2005),
		and \kbd{"Wald"} the Wald interval, which cannot be recommended.
 }

}
\details{

This function allows computation of confidence intervals for binomial group testing s>1 and for simple binomial testing if s=1.
If an actual confidence level greater or equal to that specified in the conf.level argument shall always be guaranteed, the exact method of 
Clopper-Pearson (method="CP") can be recommended for one-sided and the improved method of Blaker (method="Blaker") can be recommended for two-sided
hypotheses. If a mean confidence level close to that specified in the argument conf.level is required, but moderate violation of this level is acceptable,
the Second-Order corrected (method="SOC"), Wilson Score (method="Score") or Agresti-Coull (method="AC") might be used, where SOC has the most symmetric
coverage and Wilson and Agresti-Coull are in tendency conservative for the upper bound and proportions close to 0 and for the lower bound and proportions
close to 1. The Wald CI might be used for large number of observations n>10000 or intermediate proportions.

  }
\value{
A list containing:
  \item{conf.int }{a confidence interval for the proportion}
  \item{estimate }{the point estimator of the proportion}

}

\references{

\emph{Schaarschmidt, F. (2007).} Experimental design for one-sided confidence intervals or hypothesis tests in binomial group testing. Communications in Biometry and Crop Science 2 (1), 32-40.
 http://agrobiol.sggw.waw.pl/cbcs/

Further:

For binomial group testing and the Wald interval: 
Swallow WH, 1985: Group testing for estimating infection rates and probabilities of disease transmission. Phytopathology Vol.75, N.8, 882-889.

For application of binomial confidence intervals in Group testing
Tebbs JM & Bilder CR, 2004: Confidence interval procedures for the probability of disease transmission in multiple-vector-transfer designs. Journal of Agricultural, Biological and Environmental Statistics, Vol.9, N.1, 75-90.

For the Blaker interval:
Blaker H (2000): Confidence curves and improved exact confidence intervals for discrete distributions. The Canadian Journal of Statistics 28 (4): 783-798.

For general discussion of the different binomial confidence intervals:
Brown LD, Cai TT, DasGupta A (2001): Interval estimation for a binomial proportion. Statistical Science 2001, Vol.16, No.2, 101-133.

For the second-order corrected interval:
Cai TT(2005): One-sided confidence intervals in discrete distributions. Journal of statistical planning and inference 131, 63-88.

 }

\seealso{ \code{\link{bgtvs}} for an exact confidence interval when designs with different group sizes are used
\code{\link{bgtTest}}:  for hypothesis tests in binomial group testing 
and the functions \code{\link{nDesign}}, \code{\link{sDesign}} and \code{\link{estDesign}} for sample size calculation in binomial group testing 
binom.test(stats)
}
\examples{

# See the example in Tebbs and Bilder (2004)
# the two.sided 95-percent 
# Clopper-Pearson as default method:

bgtCI(n=24,Y=3,s=7)
bgtCI(n=24,Y=3,s=7,conf.level=0.95,
 alternative="two.sided", method="CP")

# other methods:
# Blaker CI is exact but shorter
# than Clopper-Pearson, only two.sided

bgtCI(n=24,Y=3,s=7, alternative="two.sided",
 method="Blaker")

# the asymptotic Wilson CI might even
# be shorter: 

bgtCI(n=24,Y=3,s=7, alternative="two.sided",
 method="Score")

# one-sided confidence intervals:

bgtCI(n=24,Y=3,s=7, alternative="less", method="CP")

# Wilson Score interval is less conservative 
bgtCI(n=24,Y=3,s=7, alternative="less", method="Score")

# the second-order corrected CI is even shorter
# in this situation:
bgtCI(n=24,Y=3,s=7, alternative="less", method="SOC")

}

\keyword{ htest }
\concept{binomial group testing}
\concept{confidence interval}
