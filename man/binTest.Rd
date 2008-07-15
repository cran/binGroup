\name{binTest}
\alias{binTest}

\title{ Calculate p-values fro hypothesis tests for a single binomial proportion  }
\description{
Calculates p values for hypothesis tests of a single binomial proportion.
Exact test, Score test and Wald test(not recommended) are available methods. 

}
\usage{
binTest(n, Y, p.hyp, alternative = "two.sided",
 method = "Exact")
}

\arguments{
  \item{n}{ integer, number of trials (number of individuals under observation)  }
  \item{Y}{ integer, number of successes ( number of individuals showing the trait of interest }
  \item{p.hyp}{ number between 0 and 1, specifying the hypothetical threshold proportion to test against }
  \item{alternative}{ character string defining the alternative hypothesis, either 'two.sided', 'less' or 'greater' }
  \item{method}{ character string defining the test method to be used:
can be one of "Exact" for an exact test corresponding to the Clopper-Pearson confidence interval, uses the exact binom.test(stats)
"Score" for a Score test, corresponding to the Wilson confidence interval
"Wald" for a Wald test corresponding to the Wald confidence interval, not recommended}

}

\value{
      A list containing:
  \item{p.val }{the p value of the test}
  \item{estimate }{the estimated proportion}

and the alternative and threshold proportion p.hyp specified in the function call.
}
\references{ Santner, TJ and Duffy, DE (1989): The statistical analysis of discrete data. Springer Verlag New York Berlin Heidelberg. Chapter 2.1.  }

\seealso{ binom.test(stats) for the exact test and corresponding confindence interval }
\examples{

# 200 seeds are taken from a seed lot.
# 2 are found to be defective.
# Is the proportion of defectives in 
# the lot less than 2 per cent? 
# H0: p >= 0.02  shall be rejected in favor of  HA: p < 0.02.
# The exact test shall be used for decision:

binTest(n=200, Y=2, p.hyp=0.02, alternative="less", method="Exact" )

}
\keyword{ htest }
\concept{ binomial}
