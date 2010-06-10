sim.g <- function (x = NULL, gshape = 20, gscale = 2, par,
    linkf = c("logit", "probit", "cloglog"),
    sample.size, group.size, sens = 1, spec = 1)
{
    if (is.null(x)) {
        x <- rgamma(n = sample.size, shape = gshape, scale = gscale)
        X <- cbind(1, x)
    }
    else {
        X <- cbind(1, x)
        sample.size <- nrow(X)
    }
    linkf <- match.arg(linkf)
    pijk <- switch(linkf, logit = plogis(X %*% par),
            probit = pnorm(X %*% par),
            cloglog = 1 - exp(-exp(X %*% par)))
    ind <- rbinom(n = sample.size, size = 1, prob = pijk)
    upper <- ceiling(sample.size/group.size)
    groupn <- rep(1:upper, each = group.size)[1:sample.size]
    num.g <- max(groupn)
    save.sum <- tapply(ind, groupn, sum)
    save.group <- as.vector(ifelse(save.sum > 0, 1, 0))
    save.obs <- rep(NA, num.g)
    for (i in 1:num.g)
         save.obs[i] <- ifelse(save.group[i] == 1, rbinom(1, 1, sens),
                1 - rbinom(1, 1, spec))
    gres <- rep(save.obs, each = group.size)[1:sample.size]
    grd <- data.frame(gres = gres, x = x, groupn = groupn, ind = ind)
    if (ncol(X) > 2)
        for (i in 2:ncol(X))
             colnames(grd)[i] <- paste("x", i - 1, sep="")
    grd
}
