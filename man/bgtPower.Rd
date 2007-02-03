\name{bgtPower}
\alias{bgtPower}
\alias{bgtPowerI}



\title{ Power to reject a hypothesis using confidence intervals in binomial group testing }
\description{ Closed calculation of the Power to reject a hypothesis in a binomial group testing
experiment using confidence intervals for decision. Closed calculation of bias of the point estimator
 for a given experimental design n, s and the true, unknown proportion. 
  
}
\usage{
bgtPower(n, s, delta, p.hyp, conf.level = 0.95, method = "CP", alternative = "two.sided")
}
\arguments{
  \item{n}{  integer, giving the number of groups i.e. assays i.e. observations, a vector of integers is also allowed }
  \item{s}{ integer, giving the common size of groups i.e. the number of individual units in each group, a vector of integers is also allowed  }
  \item{delta}{ absolute difference between the threshold and the true proportion which shall be detectable with specified power, a vector is also allowed  }
  \item{p.hyp}{ number between 0 and 1, specifying the threshold proportion in the hypotheses }
  \item{conf.level}{ Confidence level required for the decision on the hypotheses }
  \item{method}{ Character string, specifying the Confidence interval method (see ?bgt.CI) to be used }
  \item{alternative}{ character string, defining the alternative hypothesis, either 'two.sided', 'less' or 'greater'
where 'less' calculates the 'power of the upper confidence limit' for a true proportion p.hyp-delta,
'greater' calculates the 'power of the lower confidece limit' for a true proportion of p.hyp+delta.
'two.sided' calculates min(power(p.hyp-delta, p.hyp+delta)) for a two.sided CI, thus can result in much lower power.
Note that coverage probability and power are not necessarily symmetric for upper and lower bound of binomial CI, 
especially for Wald, Wilson Score and Agresti-Coull CI. }
}
\details{
The power of a confidence interval here is defined as the probability that a confidence interval or limit excludes the threshold parameter (p.hyp) of the null hypothesis. 
F.e. the null hypothesis H0: p>= 0.005 might be rejected, if an upper confidence limit for p does not contain p.hyp=0.005.
In case that a delta of 0.002 shall be detectable, this function calculates the probability, that an interval of a given method will
exclude p.hyp=0.005 if the true proportion = 0.003 for a given group testing design with n groups, each containing s units.
Due to discreteness, the power does not increase monotone for increasing  number of groups n or group size s, but exhibits local maxima and minima, depending on n,s, conf.level, p.hyp.
The power can be identical for different methods, depending on the particular combination of n, s, p.hyp, conf.level.

In case that vectors are specified for n,s, and (or) delta, a matrix will be constructed and power and bias are calculated for each line in this matrix. 
}
\value{
    A matrix containing the columns
  \item{n}{a vector of number of groups}
  \item{s}{a vector of group sizes}
  \item{ns}{a vector of total sample size n*s resulting from the latter}
  \item{delta}{a vector of delta}
  \item{power}{ the power to reject the given null hypothesis, with the specified method and parameters of the first 4 columns}
  \item{bias}{the bias of the estimator for the specified n, s, and the true proportion}
}
\references{
Swallow WH, 1985: Group testing for estimating infection rates and probabilities of disease transmission. Phytopathology Vol.75, N.8, 882-889.

Tebbs JM & Bilder CR, 2004: Confidence interval procedures for the probability of disease transmission in multiple-vector-transfer designs. Journal of Agricultural, Biological and Environmental Statistics, Vol.9, N.1, 75-90.
}

\author{  Frank Schaarschmidt, Bioinformatics Unit, University of Hannover, Germany  }


\seealso{ 
\code{\link{nDesign}}:  stepwise increasing n (for a fixed group size s) until a certain power is reached within a restriction of bias of the estimator
\code{\link{sDesign}}:  stepwise increasing s (for a fixed number of groups) until a certain power is reached within a restriction of bias of the estimator 
\code{\link{estDesign}}: selection of an appropriate design to achieve minimal mean square error of the estimator
 }
\examples{

# Calculate the power for the design
# in the example given in Tebbs and Bilder(2004):
# n=24 groups each containing 7 insects
# if the true proportion of virus vectors
# in the population is 0.04 (4 percent),
# the power to reject H0: p>=0.1 using an
# upper Clopper-Pearson ("CP") confidence interval
# can be calculated using the following call:

bgtPower(n=24, s=7, delta=0.06, p.hyp=0.1,
 conf.level=0.95, alternative="less", method="CP")


# c(), seq() or rep() might be used to explore development
# of power and bias for varying n, s, delta. How much can
# we decrease the number of groups (costly assays to be performed)
# by pooling the same number of 320 individuals to groups of
# increasing size without largely decreasing power ??
#

bgtPower(n=c(320,160,80,64,40,32,20,10,5),
 s=c(1,2,4,5,8,10,16,32,64),
 delta=0.01, p.hyp=0.02)

# How power develops for increasing difference
# delta of the true proportion to the threshold proportion?


bgtPower(n=50, s=10, delta=seq(from=0, to=0.01, by=0.001),
 p.hyp=0.01, method="CP")

# use a more liberal method:

bgtPower(n=50, s=10, delta=seq(from=0, to=0.01, by=0.001),
 p.hyp=0.01, method="SOC")

}

\keyword{ htest }
\concept{ sample size }
\concept{ group testing }
\concept{ binomial }
\concept{ power }
