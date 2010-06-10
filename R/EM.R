EM <-
function (Y, X, groupn, sens, spec, linkf, start = NULL, control = EM.control())
{
    z <- tapply(Y, groupn, tail, n = 1)
    K <- ncol(X)
    if (is.null(start)) {
        if (K == 1)
            cova.mean <- as.matrix(tapply(X, groupn, mean))
        else {
            temp <- by(X, groupn, mean)
            cova.mean <- do.call(rbind, temp)
        }
        beta.old <- lm.fit(cova.mean, z)$coefficients
    } 
    else beta.old <- start
    sam <- length(Y)
    vec <- 1:sam
    group.sizes <- tapply(Y, groupn, length)
    diff <- 1
    counts <- 1
    extra.loop <- FALSE
    while ((diff > control$tol || extra.loop) && counts <= control$maxit) {
        xib <- X %*% beta.old
        pijk <- switch(linkf, logit = plogis(xib),
            probit = pnorm(xib), cloglog = 1 -
                exp(-exp(xib)))
        prodp <- tapply(1 - pijk, groupn, prod)
        den <- rep((1 - spec) * prodp + sens * (1 - prodp),
            group.sizes)
        den2 <- rep(spec * prodp + (1 - sens) * (1 - prodp),
            group.sizes)
        expect <- rep(NA, times = sam)
        for (i in vec) {
            if (Y[i] == 0)
                expect[i] <- (1 - sens) * pijk[i]/den2[i]
            else expect[i] <- sens * pijk[i]/den[i]
        }
        if (!extra.loop) {
            suppress <- function(w) 
                if(any(grepl("non-integer #successes in a binomial glm", w))) 
                   invokeRestart("muffleWarning")
            mod.fit <- withCallingHandlers(glm.fit(X, expect, 
                family = binomial(link = linkf)), warning = suppress)
            diff <- max(abs((beta.old - mod.fit$coefficients)/beta.old))
            beta.old <- mod.fit$coefficients
            if (control$trace)
                cat("beta is", beta.old, "\tdiff is", diff, "\tIterations -", counts, "\n")
        }
        counts <- counts + 1
        if (!extra.loop && diff <= control$tol) {
           extra.loop <- TRUE
           count2 <- counts
           counts <- control$maxit
        }
    }    
    erf <- 2 * pijk - 1
    pt1 <- switch(linkf, logit = -exp(xib)/(1 + exp(xib))^2,
        probit = sqrt(2) * xib * exp(-xib^2/2)/(sqrt(pi) * (1 -
            erf)) - 2 * exp(-xib^2)/(pi * (1 - erf)^2), cloglog = -exp(xib))
    pt2 <- switch(linkf, logit = 0, probit = (8 * exp(-xib^2/2) *
        erf + 2 * xib * sqrt(2 * pi) * erf^2 - 2 * xib * sqrt(2 *
        pi)) * exp(-xib^2/2)/((1 + erf)^2 * pi * (1 - erf)^2),
        cloglog = -(exp(xib - exp(xib)) + exp(2 * xib - exp(xib)) -
            exp(xib))/(exp(-exp(xib)) - 1)^2)
    nm <- pt1 + expect * pt2
    sign1 <- as.vector(sign(nm))
    nn <- as.vector(sqrt(abs(nm)))
    x2 <- X * nn
    m <- (t(x2) %*% (sign1 * x2))
    b <- array(NA, c(K, K, sum(group.sizes^2)))
    p <- 1
    for (i in vec) for (j in vec[groupn == groupn[i]]) {
        wii <- ifelse(i == j, expect[i] - expect[i]^2, expect[i] *
            (pijk[j] - expect[j]))
        coe <- switch(linkf, logit = 1, probit = 8 * exp(-(xib[i]^2 +
            xib[j]^2)/2)/((1 - erf[i]^2) * (1 - erf[j]^2) * pi),
            cloglog = exp(xib[i] + xib[j])/((exp(-exp(xib[i])) -
                1) * (exp(-exp(xib[j])) - 1)))
        b[, , p] <- wii * coe * X[i, ] %*% t(X[j, ])
        p <- p + 1
    }
    m1 <- apply(b, c(1, 2), sum)
    H <- -(m + m1)
    zhat <- sens + (1 - sens - spec) * prodp
    residual <- z - zhat
    residd <- -2 * sum(z * log(zhat) + (1 - z) * log(1 - zhat))
    logL0 <- function(beta) {
        inter <- rep(beta, sam)
        pijk <- switch(linkf, logit = plogis(inter), probit = pnorm(inter),
            cloglog = 1 - exp(-exp(inter)))
        prodp <- tapply(1 - pijk, groupn, prod)
        -sum(z * log(sens + (1 - sens - spec) * prodp) +
            (1 - z) * log(1 - sens - (1 - sens - spec) * prodp))
    }
    mod.fit0 <- optim(par = .Call("logit_link", mean(z), PACKAGE = "stats"), fn = logL0,
        method = "BFGS", control = list(trace = 0, maxit = 1000))
    nulld <- 2 * mod.fit0$value
    aic <- residd + 2 * K
    if (count2 == control$maxit) count2 <- "did not converge"
    list(coefficients = beta.old, hessian = H, fitted.group.values = zhat,
        deviance = residd, df.residual = sam - K, null.deviance = nulld,
        df.null = sam - 1, aic = aic, counts = count2, residuals = residual,
        z = z)
}

