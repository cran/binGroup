print.summary.gt <-
function(x, digits=max(3, getOption("digits") - 3), 
           signif.stars = getOption("show.signif.stars"), ...)
{
    obj <- x
    cat("\nCall:\n")
    cat(paste(deparse(obj$call), sep = "\n", collapse = "\n"), 
           "\n\n", sep = "")

    cat("Deviance Residuals: \n")
    if (obj$df.residual > 5) {
           obj$deviance.resid <- quantile(obj$deviance.resid, na.rm = TRUE)
           names(obj$deviance.resid) <- c("Min", "1Q", "Median", "3Q", 
                                            "Max")
    }
    print.default(obj$deviance.resid, digits = digits, na.print = "", 
            print.gap = 2)

    cat("\nCoefficients:\n")
    coefs <- obj$coefficients
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
               na.print = "NA", ...)
    cat("\n",apply(cbind(paste(format(c("Null", 
              "Residual"), justify = "right"), "deviance:"), format(unlist(obj[c("null.deviance", 
              "deviance")]), digits = 4)," on", 
              format(unlist(obj[c("df.null", "df.residual")])), " degrees of freedom\n"), 
              1, paste, collapse = " "), sep = "")

    if (obj$method=="Vansteelandt")
    cat("AIC: ", format(obj$aic, digits = 4), 
          "\n\n", "Number of iterations in optim(): ", obj$counts, 
          "\n", sep = "")
    else 
    cat("AIC: ", format(obj$aic, digits = 4), 
          "\n\n", "Number of iterations in EM: ", obj$counts, 
          "\n", sep = "")
    cat("\n")
    invisible(obj)
}

