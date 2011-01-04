\name{gt.control}
\alias{gt.control}
\title{
Auxiliary for Controlling Group Testing Regression
}
\description{
Auxiliary function to control fitting parameters of EM algorithm used internally in \code{\link{gtreg.mp}} and \code{\link{EM.mp}} or \code{\link{gtreg}} and \code{\link{gtreg.fit}} with method = "Xie".
}
\usage{

gt.control(tol = 0.005, n.gibbs = 1000, n.burnin = 20,
 maxit = 500, trace = FALSE, time = TRUE)

}

\arguments{
  \item{tol}{
convergence criterion 
}
  \item{n.gibbs}{
the Gibbs sample size to be used in each E step for the EM algorithm (default is 1000), for matrix pooling or simple pooling with retests (Dorfman's procedure)
}
  \item{n.burnin}{
the number of samples in the burn-in period (default is 20), for matrix pooling or simple pooling with retests (Dorfman's procedure) 
}
  \item{maxit}{
maximal number of iterations in the EM algorithm
}
  \item{trace}{
logical indicating if output should be printed for each iteration, defaults to \kbd{FALSE}
}
  \item{time}{
logical indicating if the length of time for the model fitting should be printed, defaults to \kbd{TRUE}}
}

\value{A list with components named as the arguments}

\examples{
# The default settings:
gt.control()

}

\keyword{misc}
