\name{bgtTest}
\alias{bgtTest}
\alias{print.bgtTest}

\title{ p values for testing binomial proportions from binomial group testing }

\description{ Calculates p values for hypothesis tests of binomial proportions 
estimated from binomial group testing experiments against a threshold proportion
in the hypotheses. Exact test, Score test and Wald test(not recommended) are available as methods. 
Assumes equal group sizes, an assay method classifying a group as positive if at least one unit in
the group is positive, individuals units randomly assigned to the groups.
 }

\usage{
bgtTest(n, Y, s, p.hyp, alternative = "two.sided", method = "Exact")
}

\arguments{
  \item{n}{ integer, number of  groups (i.e. assays i.e. observations)  }
  \item{Y}{ integer, number of positive groups }
  \item{s}{ integer, common size of groups i.e. the number of individual units in each group}
  \item{p.hyp}{ number between 0 and 1, specifying the hypothetical threshold proportion to test against }
  \item{alternative}{ character string defining the alternative hypothesis, either 'two.sided', 'less' or 'greater' }
  \item{method}{ character string defining the test method to be used:
can be one of "Exact" for an exact test corresponding to the Clopper-Pearson confidence interval
"Score" for a Score test, corresponding to the Wilson confidence interval
"Wald" for a Wald test corresponding to the Wald confidence interval, not recommended}
}

\details{
}

\value{A list containing:
  \item{p.val }{the p value of the test}
  \item{estimate }{the estimated proportion}
  
}
\references{ 
Swallow WH, 1985: Group testing for estimating infection rates and probabilities of disease transmission. Phytopathology Vol.75, N.8, 882-889.

Blyth, C and Still, H. (1983): Binomial confidence intervals. Journal of the American Statistical Association 78, 108-116.

Santner TJ and Duffy DE (1989): The statistical analysis of discrete data. Springer Verlag New York Berlin Heidelberg.

Remund KM, Dixon DA, Wright DL, Holden LR. Statistical considerations on seed purity testing on transgenic traits. Seed Sci Res 2001; 11: 101-119.

 }
\author{  Frank Schaarschmidt, Bioinformatics Unit, University of Hannover, Germany }
\note{ }

\seealso{ \code{\link{bgtCI}} for confidence intervals in binomial group testing}
\examples{

# Assume the experiment: You perform assays on
# n=10 groups, each group is a bulk sample
# of s=100 individuals, your aim is to show the 
# 'less than 0.5 percent'('p<0.005') of the units
# of the population show a detrimental trait (positive assay).
# Your assay is senstive to show a positive result if only 1
# unit in the bulk sample of 100 units is positive.
# You observed Y=1 positive assay, 9 negative assays.

bgtTest(n=10,Y=1,s=100,alternative="less",method="Exact",p.hyp=0.005)

# Interval bounds of Clopper-Pearson
# in the example of Tebbs and Bilder(2004):

bgtTest(n=24, Y=3, s=7, alternative="two.sided",
 method="Exact", p.hyp=0.0543)

bgtTest(n=24, Y=3, s=7, alternative="two.sided",
 method="Exact", p.hyp=0.0038)

# Further methods:

bgtTest(n=24, Y=3, s=7, alternative="two.sided",
 method="Score", p.hyp=0.0516)

bgtTest(n=24, Y=3, s=7, alternative="two.sided",
 method="Wald", p.hyp=0.0401)
}

\keyword{ htest }
\concept{ binomial }
\concept{ group testing }