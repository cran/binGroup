\name{sim.mp}
\Rdversion{1.1}
\alias{sim.mp}

\title{
Simulation Function for Group Testing Data with Matrix Pooling Design
}
\description{
Simulates data in group testing form ready to be fit by \kbd{gtreg.mp}.
}
\usage{
sim.mp(gshape = 20, gscale = 2, beta.par, rown, coln,
 sens = 1, spec = 1, sens.ind = NULL, spec.ind = NULL)
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
  \item{rown}{
a vector that specifies the number of rows in each matrix, a scalar if only one matrix is simulated
}
  \item{coln}{
a vector that specifies the number of columns in each matrix, a scalar if only one matrix is simulated
}
  \item{sens}{
sensitivity of the group tests, set to be 1 by default.
}
  \item{spec}{
specificity of the group tests, set to be 1 by default.
}
  \item{sens.ind}{
sensitivity of the individual retests, set to be equal to sens if not specified otherwise.
}
  \item{spec.ind}{
specificity of the individual retests, set to be equal to spec if not specified otherwise.
}
}

\details{

\kbd{sim.mp} generates group testing data in matrix pooling form. To begin, the covariates are generated from a gamma distribution with given \kbd{gshape} and \kbd{gscale}. The individual probabilities are calculated with these covariates and the logit link using coefficients given in \kbd{beta.par}. The true individual responses are simulated next using Bernoulli distributions with these corresponding individual probabilities. The individuals are organized into (by column) one or more matrices specified by \kbd{rown} and \kbd{coln}, and the true group responses are found (i.e., if at least one response is positive, the group is positive; otherwise, the group response is negative). The row and column group responses are simulated from Bernoulli distributions using the given \kbd{sens} and \kbd{spec} values. Results of individual retests are simulated with \kbd{sens.ind} and \kbd{spec.ind} for individuals that lie on the intersection of an observed positive row and an observed positive column. 
In the case where no column (row) tests positive in a matrix, all the individuals in any observed positive rows (columns) will be assigned a simulated retest result.

}

\value{

\kbd{sim.mp} returns a list with the components
\kbd{dframe}: the data frame that is actually to be fit,
\kbd{ind}: the true individual responses presented in matrices and
\kbd{prob}: the individual probabilities.

\kbd{dframe} is a data frame with columns
\item{col.resp}{the column group response}
\item{row.resp}{the row group response}
\item{x}{the covariate}
\item{sqn}{the array number}
\item{coln}{the column group number}
\item{rown}{the row group number}
\item{retest}{the results of individual retests}

}

\author{
Boan Zhang
}

\seealso{
 \code{\link{gtreg.mp}} for the corresponding function to fit the model.
}

\examples{

# 5*6 and 4*5 matrix
set.seed(9128)
sa1a<-sim.mp(beta.par=c(-7,0.1), rown=c(5,4), coln=c(6,5),
 sens=0.95, spec=0.95)
sa1<-sa1a$dframe

}

\keyword{datagen}
\concept{group testing}
\concept{matrix pooling}
