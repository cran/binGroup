gtreg <-
function(formula, data, groupn, sens = 1, spec = 1, 
        linkf = c("logit", "probit", "cloglog"), 
        method = c("Vansteelandt", "Xie"), start = NULL, control = EM.control(...), ...) {

    call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "groupn"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    gr <- model.extract(mf, "groupn")
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm))
            names(Y) <- nm
    }
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf)
    else matrix(, NROW(Y), 0)
    linkf <- match.arg(linkf)
    fit <- switch(method <- match.arg(method),
          "Vansteelandt" = gtreg.fit(Y, X, gr, sens, spec, linkf, start),
          "Xie" = EM(Y, X, gr, sens, spec, linkf, start, control))
    fit <- c(fit, list(call = call, formula = formula, method = method, 
          link = linkf, terms = mt))
    class(fit) <- "gt"
    fit
        
}

