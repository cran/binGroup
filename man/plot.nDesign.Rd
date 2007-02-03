\name{plot.nDesign}
\alias{plot.nDesign}

\title{ Plot iteration of nDesign}
\description{
Plot results of function nDesign
}
\usage{
}

\arguments{
  \item{x}{ An object of class "nDesign" as can be created by \code{\link{nDesign}}}
  \item{\dots}{ further arguments to be passed to \code{\link{plot}} }
}

\author{ Frank Schaarschmidt, Bioinformatics Unit, University of Hannover, Germany}

\examples{

plot(nDesign(nmax=100, s=30,
 delta=0.004, p.hyp=0.005,
 alternative="less", method="Score",
 power=0.8))

}
\keyword{ hplot }
\concept{ group testing}
\concept{ sample size }