% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Array.Functions.R
\name{Informative.array.prob}
\alias{Informative.array.prob}
\title{Arrange a matrix of probabilities for informative array testing}
\usage{
Informative.array.prob(prob.vec, nr, nc, method = "sd")
}
\arguments{
\item{prob.vec}{vector of individual risk probabilities, of 
length nr*nc.}

\item{nr}{the number of rows in the array.}

\item{nc}{the number of columns in the array.}

\item{method}{character string defining the method to be used
for matrix arrangement. Options include spiral ("\kbd{sd}")
and gradient ("\kbd{gd}") arrangement. See McMahan et al. (2012)
for additional details.}
}
\value{
A matrix of probabilities arranged according to the 
specified method.
}
\description{
Arrange a vector of individual risk probabilities in a
matrix for informative array testing without master pooling.
}
\examples{
# Use the gradient arrangement method to create a matrix
#   of individual risk probabilities for a 10x10 array.
# Depending on the specified probability, alpha level, 
#   and overall group size, simulation may be necessary in 
#   order to generate the vector of individual probabilities. 
#   This is done using the p.vec.func() function and requires 
#   the user to set a seed in order to reproduce results.
# This example takes less than 1 second to run.
# Estimated running time was calculated using a 
#   computer with 16 GB of RAM and one core of an 
#   Intel i7-6500U processor.
set.seed(1107)
p.vec1 <- p.vec.func(p=0.05, alpha=2, grp.sz=100)
Informative.array.prob(prob.vec=p.vec1, nr=10, nc=10, method="gd")

# Use the spiral arrangement method to create a matrix
#   of individual risk probabilities for a 5x5 array.
# Depending on the specified probability, alpha level, 
#   and overall group size, simulation may be necessary in 
#   order to generate the vector of individual probabilities. 
#   This is done using the p.vec.func() function and requires 
#   the user to set a seed in order to reproduce results.
# This example takes less than 1 second to run.
# Estimated running time was calculated using a 
#   computer with 16 GB of RAM and one core of an 
#   Intel i7-6500U processor.
set.seed(8791)
p.vec2 <- p.vec.func(p=0.02, alpha=0.5, grp.sz=25)
Informative.array.prob(prob.vec=p.vec2, nr=5, nc=5, method="sd")
}
\references{
\insertRef{McMahan2012b}{binGroup}
}
\seealso{
\code{\link{p.vec.func}} for generating a vector of individual
risk probabilities from an overall probability of disease and
\code{\link{Array.Measures}} for calculating operating characteristics
for informative array testing without master pooling.

\url{http://chrisbilder.com/grouptesting/}

Other Individual risk probability functions: 
\code{\link{beta.dist}()},
\code{\link{p.vec.func}()}
}
\author{
This function was originally written by Christopher S. McMahan 
for McMahan et al. (2012). The function was obtained from 
\url{http://chrisbilder.com/grouptesting/}.
}
\concept{Individual risk probability functions}
