\name{plot.sDesign}
\alias{plot.sDesign}

\title{ Plot results of sDesign }
\description{
Plot function for objects of class "sDesign"
}
\usage{
\method{plot}{sDesign}(x, \dots)
}

\arguments{
  \item{x}{ An object of class "sDesign", as can be created by \code{\link{sDesign}} }
  \item{\dots}{ further arguments to be passed to \code{\link{plot}} }
}

\value{ A plot.}


\seealso{
\code{\link{bgtPower}}: calculation of power and bias depending on n, s, delta, p.hyp, conf.level, method
\code{\link{nDesign}}:  function to iteratively optimize sample size(number of groups)n for a given group size s
\code{\link{estDesign}}:  function to choose group size s according to the minimal mse of the estimator, as given in Swallow (1985)
}

\examples{

plot( sDesign(smax=100, n=30,
 delta=0.004, p.hyp=0.005,
 alternative="less", method="Score",
 power=0.8))

}
\keyword{ hplot }
\concept{ group testing}
\concept{ sample size }
\concept{ power}

