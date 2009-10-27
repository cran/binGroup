summary.gt <-
function(object, ...) {
    
    coef.p<-object$coefficients
    cov.mat<-solve(object$hessian)
    dimnames(cov.mat)<-list(names(coef.p),names(coef.p))
    var.cf<-diag(cov.mat)
    s.err<-sqrt(var.cf)
    zvalue<-coef.p/s.err
    dn <- c("Estimate", "Std. Error")
    pvalue <- 2 * pnorm(-abs(zvalue))
    coef.table <- cbind(coef.p, s.err, zvalue, pvalue)
    dimnames(coef.table) <- list(names(coef.p), c(dn, "z value", "Pr(>|z|)"))
    keep <- match(c("call", "link", "aic", "deviance",
           "df.residual", "null.deviance", "df.null",
           "counts", "method"), names(object), 0)

    ans<-c(object[keep],list(coefficients = coef.table, deviance.resid = residuals(object, type="deviance"),
           cov.mat =cov.mat))
    class(ans)<-"summary.gt"
    return(ans)
    
}

