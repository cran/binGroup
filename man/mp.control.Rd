\name{mp/EM.control}
\alias{mp.control}
\alias{EM.control}
\title{
Control arguments for fitting parameters functions gtreg and gtreg.mp
}
\description{
Auxilliary function to control the fitting parameters used internally in \code{\link{gtreg}} and \code{\link{EM}} and for the matrix pooling functions \code{\link{gtreg.mp}} and \code{\link{EM.mp}}.
}
\usage{

EM.control(tol = 0.001, maxit = 1000, trace = FALSE)

mp.control(tol = 0.005, n.gibbs = 1000, n.burnin = 20,
 maxit = 500, trace = FALSE, time = TRUE)

}

\arguments{
  \item{tol}{
convergence criterion 
}
  \item{n.gibbs}{
the Gibbs sample size to be used in each E step for the EM algorithm (default is 1000).
}
  \item{n.burnin}{
the number of samples in the burn-in period (default is 20).
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

\value{
A list with components named as the arguments of \kbd{EM.control} and \kbd{mp.control}, respectively.
}

\examples{
# The default settings:
EM.control()
mp.control()
}

\keyword{misc}
