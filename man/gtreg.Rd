\name{gtreg}
\Rdversion{1.1}
\alias{gtreg}
\alias{gtreg.fit}
\alias{EM}

\title{
Fitting Group Testing Models
}

\description{
\kbd{gtreg} is a function to fit the group testing regression model specified through a symbolic description of the linear predictor and descriptions of the group testing setting.
}

\usage{
gtreg(formula, data, groupn, sens = 1,
 spec = 1, linkf = c("logit", "probit", "cloglog"),
 method = c("Vansteelandt", "Xie"), ...)

gtreg.fit(resp, cova, groupn, sens, spec,
 linkf, optim.meth = "Nelder-Mead")

EM(resp, cova, groupn, sens, spec, linkf,
 EM.maxiter = 1000)
}


\arguments{
  \item{formula}{
an object of class \kbd{"formula"} (or one that can be coerced to that class): a symbolic description of 
the model to be fitted. The details of model specification are given under 'Details'.

}
  \item{data}{
an optional data frame, list or environment (or object coercible by \kbd{as.data.frame} to a data 
frame) containing the variables in the model. If not found in data, the variables are taken from 
\kbd{environment(formula)}, typically the environment from which \kbd{gtreg} is called.

}
  \item{groupn}{
a vector, list or data frame of the group numbers that designates individuals to groups.
}
  \item{sens}{
sensitivity of the test, set to be 1 by default.
}
  \item{spec}{
specificity of the test, set to be 1 by default.
}
  \item{linkf}{
a character string specifying one of the three link functions for a binomial model: \kbd{"logit"} 
(default) or \kbd{"probit"} or \kbd{"cloglog"}.

}
  \item{method}{
The method to fit the model, must be one of \kbd{"Vansteelandt"} (default) or \kbd{"Xie"}. The option 
\kbd{"Vansteelandt"} finds estimates by directly maximizing the likelihood function based on the group responses while the option \kbd{"Xie"} uses the EM algorithm to maximize the likelihood function in terms of the unobserved individual responses.
}

\item{resp}{
For \kbd{gtreg.fit} and \kbd{EM}: the vector of the group response variable
}

\item{cova}{
For \kbd{gtreg.fit} and \kbd{EM}: the design matrix of the covariates
}

\item{optim.meth}{
For \kbd{gtreg.fit}: the method in optim to fit the full model (default is \kbd{"Nelder-Mead"})
}

\item{EM.maxiter}{
For \kbd{EM}: The maximal number of iterations in the EM algorithm if \kbd{"Xie"} is chosen.
}

\item{\dots}{
In \kbd{gtreg}: arguments to be passed to \kbd{gtreg.fit} or \kbd{EM}.
}

}
\details{

A typical predictor has the form \kbd{groupresp ~ covariates} where response is the (numeric) group response vector and \kbd{covariates} is a series of terms which specifies a linear predictor for individual responses. Note that it is actually the unobserved individual responses, not the observed group responses, which are modeled by the covariates here. In \kbd{groupresp}, a 0 denotes a negative response and a 1 denotes a positive response, where the probability of an individual positive response is being modeled directly.  A terms specification of the form \kbd{first + second} indicates all the terms in \kbd{first} together with all the terms in \kbd{second} with duplicates removed. The \kbd{terms} in the formula will be re-ordered so that main effects come first, followed by the interactions, all second-order, all third-order and so on; to avoid this pass a terms object as the formula.

A specification of the form \kbd{first:second} indicates the set of terms obtained by taking the interactions of all terms in \kbd{first} with all terms in \kbd{second}. The specification \kbd{first*second} indicates the cross of \kbd{first} and \kbd{second}. This is the same as \kbd{first + second + first:second}. 

Two workhorse functions \kbd{gtreg.fit} and \kbd{EM}, which corresponds to Vansteelandt's and Xie's method respectively, are called by \kbd{gtreg} to carry out the model fitting. The \kbd{gtreg.fit} function uses the \kbd{optim} function to maximize the likelihood function written in terms of the group responses. The \kbd{EM} function applies Xie's EM algorithm to the likelihood function written in terms of the unobserved individual responses; the function uses \kbd{optim} to maximize the likelihood function in each M step. The EM algorithm usually converges more slowly than the method of Vansteelandt. There could be slight differences in the estimates between the two methods due to different convergence criteria.

Note the data used here should be in the form of simple pooling - meaning that each individual appears in exactly one pool. Also, no individual tests or retests are used in the model fitting. 

For the background on the use of \kbd{optim}, see \kbd{help(optim)}.

}
\value{
\kbd{gtreg} returns an object of class \kbd{"gt"}. See later in this section. 
The function \kbd{summary} (i.e., \code{\link{summary.gt}}) can be used to obtain or print a summary of the results. 
The group testing functions \kbd{predict} (i.e., \code{\link{predict.gt}}) and \kbd{residuals} (i.e., \code{\link{residuals.gt}}) can be used to extract various useful features of the value returned by \kbd{gtreg}. 
An object of class \kbd{"gt"} is a list containing at least the following components: 

 \item{ coefficients }{ a named vector of coefficients }
 \item{ hessian }{ estimated Hessian matrix of the negative log likelihood function, as an estimator of the information matrix }
 \item{ residuals }{ the response residuals, difference of the observed group responses and the fitted group responses. }
 \item{ fitted.group.values }{ the fitted mean values of group responses. }
 \item{ deviance }{ the deviance between the fitted model and the saturated model. }
 \item{ aic }{ Akaike's An Information Criterion, minus twice the maximized log-likelihood plus twice the number of coefficients }
 \item{ null.deviance }{ The deviance for the null model, comparable with \kbd{deviance}. The null model will include only the intercept if there is one in the model. }
 \item{ counts }{ For Vansteelandt's method: the number of iterations in \kbd{optim}; For Xie's method: the number of iterations in the EM algorithm.}
 \item{ df.residual }{ the residual degrees of freedom. }
 \item{ df.null }{ the residual degrees of freedom for the null model. }
 \item{ z }{ the vector of group responses. }
 \item{ call }{ the matched call. }
 \item{ formula }{ the formula supplied. }
 \item{ terms }{ the terms object used. }
 \item{ method }{ the method (\kbd{"Vansteelandt"} or \kbd{"Xie"}) used to fit the model. }
 \item{ link }{ the link function used in the model. }
}
\references{

\emph{Xie, M. (2001)}, Regression analysis of group testing samples, \emph{Statistics in Medicine, 20, 1957-1969.}

\emph{Vansteelandt, S., Goetghebeur, E., and Verstraeten, T. (2000)}, Regression models for disease prevalence with diagnostic tests on pools of serum samples, \emph{Biometrics, 56, 1126-1133.}

}

\author{
Boan Zhang
}


\seealso{

\code{\link{summary.gt}}, \code{\link{predict.gt}} and \code{\link{residuals.gt}} for \kbd{gt} methods.
\code{\link{gtreg.mp}} for the group testing regression model in the matrix pooling setting.


}
\examples{

data(hivsurv)

fit1 <- gtreg(formula = groupres ~ AGE + EDUC., data = hivsurv,
           groupn = gnum, sens = 0.9, spec = 0.9, method = "Xie")

## --- Continuing the Example from  '?sim.g':

set.seed(1125)
gt.data<-sim.g(beta.par=c(-9,0.1), number.sample=279, group.size=9)
fit2 <- gtreg(formula = groupres ~ x, data = gt.data, 
           groupn = gnum, linkf = "probit")

}

\keyword{models}
\keyword{regression}
\concept{group testing}