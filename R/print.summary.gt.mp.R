print.summary.gt.mp <-
function(x, digits=max(3, getOption("digits") - 3), 
           signif.stars = getOption("show.signif.stars"), ...)
{
    obj <- x
    cat("\nCall:\n")
    cat(paste(deparse(obj$call), sep = "\n", collapse = "\n"), 
           "\n\n", sep = "")

    cat("\nCoefficients:\n")
    coefs <- obj$coefficients
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
               na.print = "NA", ...)
    cat("\nNull Degrees of Freedom:", obj$df.null, "\nResidual Degrees of Freedom:", 
        obj$df.residual, "\n")

    cat("\nNumber of Gibbs samples generated in each E step: ", obj$Gibbs.sample.size, 
          "\n", "Number of iterations in EM algorithm: ", obj$counts, 
          "\n", sep = "")
    cat("\n")
    invisible(obj)
}

