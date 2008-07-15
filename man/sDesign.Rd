\name{sDesign}
\alias{sDesign}


\title{ Find a group size for which power achieves a specified level  }
\description{
Increasing of group size s for a fixed number of groups n in a binomial group testing design,
until a pre-specified power is achieved. Control of bias during iteration. A hypothetical
threshold proportion p.hyp and the absolute difference delta to be detected have to be specified. 
 }
\usage{
sDesign(n, smax, delta, p.hyp, conf.level = 0.95,
 power = 0.8, alternative = "two.sided", method = "CP", biasrest = 0.05)
}

\arguments{
  \item{n}{ integer, fixed sample size (number of assays) }
  \item{smax}{ integer, maximal group size allowed in planning of the design }
  \item{delta}{ absolute difference between the threshold and the true proportion which shall be detectable with specified power }
  \item{p.hyp}{ threshold proportion to test against in the hypothesis, specify as a value between 0 and 1 }
  \item{conf.level}{ Confidence level of the decision, default is 0.95 }
  \item{power}{ Level of power to be achieved to be specified as a probability between 0 and 1 }
  \item{alternative}{ character string, defining the alternative hypothesis, either 'two.sided', 'less' or 'greater'
where 'less' calculates the 'power of the upper confidence limit' for a true proportion p.hyp-delta,
'greater' calculates the 'power of the lower confidece limit' for a true proportion of p.hyp+delta.
'two.sided' calculates min(power(p.hyp-delta, p.hyp+delta)) for a two.sided CI, thus can result in much lower power.
Note that coverage probability and power are not necessarily symmetric for upper and lower bound of binomial CI, 
especially for Wald, Wilson Score and Agresti-Coull CI. }

  \item{method}{ character string specifying the CImethod to be used for evaluation, see argument method in bgtCI }
  \item{biasrest}{ value between 0 and 1 specifying the absolute bias maximally allowed }
}
\details{
The power of a confidence interval here is defined as the probability that a confidence interval or limit excludes the threshold parameter (p.hyp) of the hypothesis. 
This function increases size of groups (number of units in a bulk sample) until a pre-specified power is reached.
Since the power does not increase monotone with increasing s for binomial proportions but oscillates between local maxima and minima, 
the simple iteration given here will generally result in selecting those s, for which the given CI method shows a local minimum of 
coverage if the null hypothesis is true.
Since the positive bias of the estimator in group testing increases with increasing group size, it is checked whether bias is smaller than a 
pre-specified value (bias.rest). If bias violates this restriction for a given combination n, s, delta, s will not be further increased the
actual power of the last acceptable group size s is returned.

 
}
\value{ A list containing:

  \item{sout }{the group size necessary to met the conditions}
  \item{powerout }{the exact power for the specified parameters and the group size}
  \item{biasout }{the bias for the specified parameters and the iterated group size}
  number of values specified in the function call or produced in the iteration, which are only necessary to apply
the function plot() on objects of class 'sDesign' 

}
\references{ 

\emph{Schaarschmidt, F. (2007).} Experimental design for one-sided confidence intervals or hypothesis tests in binomial group testing. Communications in Biometry and Crop Science 2 (1), 32-40.
 http://agrobiol.sggw.waw.pl/cbcs/

For closed calculation of power see the closed calculation of coverage probability in: 
Tebbs JM & Bilder CR, 2004: Confidence interval procedures for the probability of disease transmission in multiple-vector-transfer designs. Journal of Agricultural, Biological and Environmental Statistics, Vol.9, N.1, 75-90.

For binomial group testing, problems of experimental design and calculation of bias see: 
Swallow WH, 1985: Group testing for estimating infection rates and probabilities of disease transmission. Phytopathology Vol.75, N.8, 882-889.


 }

\seealso{
\code{\link{plot.sDesign}} can be used to plot th iteration of this function

\code{\link{bgtPower}}: calculation of power and bias depending on n, s, delta, p.hyp, conf.level, method
\code{\link{nDesign}}:  function to iteratively optimize sample size(number of groups)n for a given group size s
\code{\link{estDesign}}:  function to choose group size s according to the minimal mse of the estimator, as given in Swallow (1985)

}

\examples{

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

sDesign(smax=100, n=30, delta=0.004, p.hyp=0.005,
 alternative="less", method="Score", power=0.8)

}
\keyword{ htest }
\concept{ group testing}
\concept{ sample size }
\concept{ power}