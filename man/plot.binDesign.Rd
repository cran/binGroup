\name{plot.binDesign}
\alias{plot.binDesign}

\title{Plot Results of binDesign}

\description{
Plot function to visualize the power curve.
}
\usage{
\method{plot}{binDesign}(x, \dots)
}
\arguments{
  \item{x}{ an object of class "binDesign, as can be created by \code{link{binDesign}}"}
  \item{\dots}{plot parameters as described in \code{\link{par}}}
}

\value{A plot. }

\seealso{ \code{\link{binPower}} for calculation of power }
\examples{

# Find a sample size for which the power to reject
# H0: p >= 0.1 in favor of HA: p < 0.1 
# is at least 0.9 (90 percent) in case that 
# the true proportion is 0.04 (i.e. an absolute delta 
# of 0.06 to the threshold proportion p.hyp=0.1)
# The exact one sided Clopper-Pearson CI shall be used
#  with default confidence level = 0.95.


sasi<-binDesign( nmax=200, delta=0.06, p.hyp=0.1,
 alternative="less", method="CP", power=0.9)

sasi

# Plot the result
 
# plot(sasi)


# for larger sample sizes this can be very time consuming.
# Better to use only a smaller range of n then:

sasi<-binDesign( nmax=c(200,300), delta=0.03, p.hyp=0.1,
 alternative="less", method="CP", power=0.9)

plot(sasi)


}
\keyword{ hplot }