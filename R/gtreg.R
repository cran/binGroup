gtreg <-
function(formula, data, groupn, sens=1, spec=1, linkf=c("logit","probit","cloglog"), 
        method=c("Vansteelandt","Xie"), ...) {

    call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "groupn"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    last<-dim(mf)[2]
    gr<-mf[,last]
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
    fit<-switch(method<-match.arg(method),
          "Vansteelandt"=gtreg.fit(Y, X, gr, sens, spec, linkf, ...),
          "Xie"=EM(Y, X, gr, sens, spec, linkf, ...))
    fit <- c(fit, list(call = call, formula = formula, method = method, link=linkf, terms=mt))
    class(fit)<-"gt"
    fit
        
}

