\name{binGroup-package}
\alias{binGroup}
\docType{package}
\title{
binGroup
}
\description{
 In this package, confidence intervals for the estimation of one binomial proportion from binomial
 group testing are provided. Additional to estimation, there are a number of functions for experimental
 design in the one-sample problem, based on the bias of the point estimate, the power of a test or the expected
 width of confidence intervals. 
}


\author{

 Frank Schaarschmidt

Maintainer: Frank Schaarschmidt <schaarschmidt@biostat.uni-hannover.de>

}
\references{

Key references are:

\emph{Hepworth G (1996).} Exact confidence intervals for proportions estimated by group testing. Biometrics 52, 1134-1146.

\emph{Schaarschmidt, F. (2007).} Experimental design for one-sided confidence intervals or hypothesis tests in binomial group testing. Communications in Biometry and Crop Science 2 (1), 32-40.
 http://agrobiol.sggw.waw.pl/cbcs/

\emph{Swallow WH (1985).} Group testing for estimating infection rates and probabilities of disease transmission. Phytopathology Vol.75, N.8, 882-889.

\emph{Tebbs JM & Bilder CR (2004).} Confidence interval procedures for the probability of disease transmission in multiple-vector-transfer designs.
 Journal of Agricultural, Biological and Environmental Statistics, Vol.9, N.1, 75-90.

For further details, see:

Schaarschmidt, F. (2005). Group testing - design and analysis. Thesis Fachbereich Gartenbau, Universitaet Hannover.
 http://www.biostat.uni-hannover.de/research/thesis

Blaker H (2000). Confidence curves and improved exact confidence intervals for discrete distributions. The Canadian Journal of Statistics 28 (4): 783-798.

Brown LD, Cai TT, DasGupta A (2001). Interval estimation for a binomial proportion. Statistical Science 2001, Vol.16, No.2, 101-128.

Cai, TT (2005). One-sided confidence intervals in discrete distributions. Journal of Statistical Planning and Inference 131: 63-88.
}

\keyword{ package }

\examples{

# # # 1)

# Confidence intervals for designs with equal group size (pool size):
# ?bgtCI

# See the example in Tebbs and Bilder (2004)
# the two.sided 95-percent 
# Clopper-Pearson as default method:

bgtCI(n=24,Y=3,s=7)
bgtCI(n=24,Y=3,s=7,conf.level=0.95,
 alternative="two.sided", method="CP")

# # # 2)

# Confidence intervals for designs with unequal group size (pool size):
# ?bgtvs

# The examples of Hepworth (1996), table 5 are:

 bgtvs(n=c(2,3), s=c(5,2), Y=c(0,0))

 bgtvs(n=c(2,3), s=c(5,2), Y=c(0,1))

# ...

# # # 3)

# For experimental design based on the bias of the point estimate,
# ?estDesign

### Compare Table 1 in Swallow(1985), p.885:

estDesign(n=10, smax=100, p.tr=0.001)

estDesign(n=10, smax=100, p.tr=0.01)


# # # 4)

# For experimental design based on the power in a hypothesis test,
# ?nDesign
# ?sDesign

## We want to show that a proportion is smaller
## 0.005 (i.e. 0.5 per cent) with a power 
## of 0.80 (i.e. 80 per cent) if the unknown proportion
## in the population is 0.003 (i.e. 0.3 per cent),
## thus we want to detect a delta of 0.002.
## The Clopper Pearson CI shall be used. 
## The maximal group size because of limited
## sensitivity of assay might be s=20 and we
## can only afford to perform maximally 100 assays:

nDesign(nmax=100, s=20, delta=0.002, p.hyp=0.005,
 alternative="less", method="CP", power=0.8)



## We want to show that a proportion is smaller
## than 0.005 (i.e. 0.5%) with a 
## power of 0.80 (i.e. 80%) if the unknown proportion
## in the population is 0.003 (i.e. 0.3%), thus a 
## delta = 0.002 shall be detected. A 95-per-cent 
## Clopper-Pearson CI (corresponding to an exact test)
## shall be used. The maximal number of groups might 
## be 30 where the assay sensitivity is not limited 
## until groupsize s = 100.

sDesign(smax=100, n=30, delta=0.002, p.hyp=0.005,
 alternative="less", method="CP", power=0.8)

## We might accept to detect delta=0.004,
## i.e. we accept to reject H0: p>=0.005 with 
## power 80 per cent when the true proportion is 0.001:

sDesign(smax=100, n=30, delta=0.004, p.hyp=0.005,
 alternative="less", method="CP", power=0.8)

# # # 

# For experimental design based on the expected width of confidence intervals,
# ?bgtWidth

# There is a minimal expected CI length, if 
# group size s is increased (fixed other parameters)
# the corresponding group size might be chosen:

bgtWidth(n=20, s=seq(from=1, to=200, by=10),
 p=0.01, alternative="less", method="CP" )

}
