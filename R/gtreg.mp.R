gtreg.mp <-
function (formula, data, coln, rown, arrayn, retest = NULL,
    sens = 1, spec = 1, linkf = c("logit", "probit", "cloglog"), 
    sens.ind = NULL, spec.ind = NULL, start = NULL, control = mp.control(...), ...)
{
    call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "coln", "rown",
        "arrayn"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    arrayn <- model.extract(mf, "arrayn")
    rown <- model.extract(mf, "rown")
    coln <- model.extract(mf, "coln")
    if (!is.na(pos <- match(deparse(substitute(retest)), names(data))))
        retest <- data[, pos]
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
    fit <- EM.mp(Y[, 1], Y[, 2], X, coln, rown, arrayn, retest,
        sens, spec, linkf, sens.ind, spec.ind, start, control)
    fit <- c(fit, list(call = call, formula = formula, link = linkf,
        terms = mt))
    class(fit) <- c("gt.mp", "gt")
    fit
}

