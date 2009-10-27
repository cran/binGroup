\name{sim.g}
\Rdversion{1.1}
\alias{sim.g}

\title{
Simulation Function for Group Testing Data
}
\description{
Simulates data in group testing form ready to be fit by \kbd{gtreg}.
}
\usage{
sim.g(gshape = 20, gscale = 2, beta.par, number.sample, 
group.size, sens = 1, spec = 1)
}

\arguments{
  \item{gshape}{
shape parameter of gamma distribution, must be non-negative, set to be 20 by default
}
  \item{gscale}{
scale parameter of gamma distribution, must be strictly positive, set to be 2 by default
}
  \item{beta.par}{
the true coefficients in the linear predictor 
}
  \item{number.sample}{
sample size of simulated data
}
  \item{group.size}{
group size in pooling individual samples
}
  \item{sens}{
sensitivity of the group tests, set to be 1 by default.
}
  \item{spec}{
specificity of the group tests, set to be 1 by default.
}
}
\details{
\kbd{sim.g} generates group testing data in simple pooling form. To begin, the covariates are generated from a gamma distribution with given \kbd{gshape} and \kbd{gscale}. The individual probabilities are calculated with these covariates and the logit link using coefficients given in \kbd{beta.par}. The true individual responses are simulated next using Bernoulli distributions with these corresponding individual probabilities. The individual responses and their corresponding covariates are ordered (by covariate) and pooled into groups by the specified \kbd{group.size}. The true group responses are found from the individual responses within the groups (i.e., if at least one response is positive, the group is positive; otherwise, the group response is negative). Finally, the observed group responses are simulated from Bernoulli distributions using the given \kbd{sens} and \kbd{spec} values. Note that with the given group size, the last group may have less individuals.

}
\value{
\kbd{sim.g} returns a data frame with the following columns:
\item{groupres}{the group response}
\item{x}{the covariate}
\item{gnum}{the group number}
\item{ind}{the actual individual response}
}

\author{
Boan Zhang
}

\seealso{
\code{\link{gtreg}}, \code{\link{gtreg.mp}}
}
\examples{

set.seed(1125)
gt.data<-sim.g(beta.par=c(-9,0.1), number.sample=279, group.size=9)

}

\keyword{datagen}
\concept{group testing}