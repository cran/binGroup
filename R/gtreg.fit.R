gtreg.fit <-
function (Y, X, groupn, sens, spec, linkf, start = NULL)
{
    z <- tapply(Y, groupn, tail, n = 1)
    K <- ncol(X)
    sam <- length(Y)
    if (is.null(start)) {
        if (K == 1) {
            cova.mean <- as.matrix(tapply(X, groupn, mean)) 
            optim.meth <- "BFGS"
        } 
        else {
            temp <- by(X, groupn, mean)
            cova.mean <- do.call(rbind, temp)
            optim.meth <- "Nelder-Mead"
        }
        beta.group <- glm.fit(cova.mean, as.vector(z), 
            family = binomial(link = linkf))$coefficients
    } 
    else {
        beta.group <- start
        names(beta.group) <- dimnames(X)[[2]]
        optim.meth <- ifelse(K == 1, "BFGS", "Nelder-Mead")
    }   
    logL <- function(beta) {
        pijk <- switch(linkf, logit = plogis(X %*% beta), probit = pnorm(X %*%
            beta), cloglog = 1 - exp(-exp(X %*% beta)))
        prodp <- tapply(1 - pijk, groupn, prod)
        -sum(z * log(sens + (1 - sens - spec) * prodp) +
            (1 - z) * log(1 - sens - (1 - sens - spec) * prodp))
    }
    mod.fit <- optim(par = beta.group, fn = logL, method = optim.meth,
        control = list(trace = 0, maxit = 1000), hessian = TRUE)
    if (det(mod.fit$hessian) == 0)
        mod.fit <- optim(par = beta.group, fn = logL, method = "SANN", hessian = TRUE)
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
    residd <- 2 * mod.fit$value
    xib <- X %*% mod.fit$par
    pijk <- switch(linkf, logit = plogis(xib), probit = pnorm(xib),
        cloglog = 1 - exp(-exp(xib)))
    prodp <- tapply(1 - pijk, groupn, prod)
    zhat <- sens + (1 - sens - spec) * prodp
    residual <- z - zhat
    aic <- residd + 2 * K
    if (mod.fit$convergence == 0)
        counts <- mod.fit$counts[[1]]
    else counts <- "did not converge"
    list(coefficients = mod.fit$par, hessian = mod.fit$hessian,
        fitted.group.values = zhat, deviance = residd, df.residual = sam - K,
        null.deviance = nulld, df.null = sam - 1, aic = aic, counts = counts,
        residuals = residual, z = z)
}

