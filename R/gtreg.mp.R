gtreg.mp <-
function(formula, data, col.groupn, row.groupn, arrayn, retest=NULL, sens=1, spec=1, 
                 linkf=c("logit","probit","cloglog"), ...) {

    call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "col.groupn", "row.groupn", "arrayn"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    last<-dim(mf)[2]
    arrayn<-mf[,last]; rown<-mf[,last-1]; coln<-mf[,last-2]
    if (!is.na(pos<-match(deparse(substitute(retest)), names(data)))) retest<-data[,pos]
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm))
            names(Y) <- nm
    }
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf)  #create design matrix
    else matrix(, NROW(Y), 0)
    linkf<-match.arg(linkf)
    fit <- EM.mp(Y[,1], Y[,2], X, coln, rown, arrayn, retest, sens, spec, linkf, ...)
    fit <- c(fit, list(call = call, formula = formula, link=linkf, terms=mt))
    class(fit)<-c("gt.mp", "gt")
    fit
        
}

