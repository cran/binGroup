\name{bgtvs}
\alias{bgtvs}
\title{ Binomial group testing with variable group sizes }
\description{
Calculates confidence intervals for binomial group testing
if groups of different size are evaluated
}
\usage{
bgtvs(n, s, Y, conf.level = 0.95, alternative = "two.sided",
 maxiter = 100)
}
\arguments{
  \item{n}{ numeric vector, specifying the number of groups of the different sizes }
  \item{s}{ numeric vector, specifying the group sizes, must be of same length as \code{n} }
  \item{Y}{ numeric vector, specifying the number of positive groups among \code{n}}
  \item{conf.level}{ a single numeric value, the confidence level of the interval }
  \item{alternative}{ a character string, with options "two.sided", "less", "greater" }
  \item{maxiter}{ maximal number steps in iteration of confidence limits  }
}
\details{
 Hepworth (1996) describes methods for constructing confidence intervals in binomial group testing,
if groups of different size are used. Currently, only the exact method (Hepworth, 1996, equation5, Table.5)
is implemented. Note, that the algorithm becomes computationally very extensive if the number of different
groups becomes larger than 3. See Hepworth (1996) for application of the method

}
\value{A list containing
  \item{conf.int }{a numerioc vector, the lower and upper limits of the confidence interval}
  \item{estimate }{the point estimate}
  moreover, some of the input arguments.
}
\references{ Hepworth, G (1996): Exact confidence intervals for proportions estimated by group testing. Biometrics 52, 1134-1146. }


\examples{
# Consider a very simple example,
# given in Hepworth (1996), table 5:
# 2 groups each containing 5 units,
# and 3 groups, each containing 2 units
# evaluated for containing at least one 
# positive unit:

# For the result Y=1 positive group among the groups
# of size 5, and Y=2 positive groups among the groups
# of size 2 we yield the interval:

bgtvs(n=c(2,3), s=c(5,2), Y=c(1,2)) 

###############################################

# The examples of Hepworth (1996), table 5 are:


 bgtvs(n=c(2,3), s=c(5,2), Y=c(0,0))

 bgtvs(n=c(2,3), s=c(5,2), Y=c(0,1))

 bgtvs(n=c(2,3), s=c(5,2), Y=c(0,2))

 bgtvs(n=c(2,3), s=c(5,2), Y=c(0,3))

 bgtvs(n=c(2,3), s=c(5,2), Y=c(1,0))

 bgtvs(n=c(2,3), s=c(5,2), Y=c(1,1))

 bgtvs(n=c(2,3), s=c(5,2), Y=c(1,2))

 bgtvs(n=c(2,3), s=c(5,2), Y=c(1,3))

 bgtvs(n=c(2,3), s=c(5,2), Y=c(2,0))

 bgtvs(n=c(2,3), s=c(5,2), Y=c(2,1))

 bgtvs(n=c(2,3), s=c(5,2), Y=c(2,2))

 bgtvs(n=c(2,3), s=c(5,2), Y=c(2,3))

}

\keyword{ htest }
\concept{group testing}
\concept{variable group size}
