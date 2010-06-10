sim.mp <-
function (x = NULL, gshape = 20, gscale = 2, par,
    linkf = c("logit", "probit", "cloglog"),
    n.row, n.col, sens = 1, spec = 1, sens.ind = NULL, spec.ind = NULL)
{
    if (is.null(sens.ind))
        sens.ind <- sens
    if (is.null(spec.ind))
        spec.ind <- spec
    if (length(n.row) != length(n.col))
        stop("vector n.row and n.col must have the same length")
    linkf <- match.arg(linkf)
    if (is.null(x)) {
        sample.size <- sum(n.col * n.row)
        x <- rgamma(n = sample.size, shape = gshape, scale = gscale)
        X <- cbind(1, x)
    }
    else {
        X <- cbind(1, x)
        sample.size <- nrow(X)
        if (sum(n.col * n.row) != sample.size)
            stop("n.row and n.col not consistent with the sample size")
    }
    len <- length(n.row)
    pijk <- switch(linkf, logit = plogis(X %*% par),
            probit = pnorm(X %*% par),
            cloglog = 1 - exp(-exp(X %*% par)))
    ind <- rbinom(n = sample.size, size = 1, prob = pijk)
    individual <- col.groupn <- row.groupn <- numeric(0)
    rowr <- colr <- numeric(0)
    ret <- rep(NA, sample.size)
    for (i in 1:len) {
        if (i > 1)
            index <- seq(max(index) + 1, length = (n.row * n.col)[i])
        else index <- 1:(n.row * n.col)[1]
        indm <- matrix(ind[index], nrow = n.row[i])
        col.resp <- apply(indm, MARGIN = 2, FUN = sum)
        col.resp <- ifelse(col.resp > 0, 1, 0)
        col.err <- rep(NA, n.col[i])
        for (j in 1:n.col[i])
             col.err[j] <- ifelse(col.resp[j] == 1, rbinom(1, 1, sens),
                                 1 - rbinom(1, 1, spec))
        row.resp <- apply(indm, MARGIN = 1, FUN = sum)
        row.resp <- ifelse(row.resp > 0, 1, 0)
        row.err <- rep(NA, n.row[i])
        for (j in 1:n.row[i])
             row.err[j] <- ifelse(row.resp[j] == 1, rbinom(1, 1, sens),
                                 1 - rbinom(1, 1, spec))
        temp.c <- rep(1:n.col[i], each = n.row[i])
        col.groupn <- c(col.groupn, temp.c)
        temp.r <- rep(1:n.row[i], n.col[i])
        row.groupn <- c(row.groupn, temp.r)
        temp2.c <- rep(col.err, each = n.row[i])
        colr <- c(colr, temp2.c)
        temp2.r <- rep(row.err, n.col[i])
        rowr <- c(rowr, temp2.r)
        if (all(row.err == 0)) {
            for (j in index) {
                 if (colr[j] == 1)
                     ret[j] <- ifelse(ind[j] == 1, rbinom(1,
                        1, sens.ind), 1 - rbinom(1, 1, spec.ind))
            }
        }
        else {
            if (all(col.err == 0)) {
                for (j in index) {
                     if (rowr[j] == 1)
                         ret[j] <- ifelse(ind[j] == 1, rbinom(1,
                            1, sens.ind), 1 - rbinom(1, 1, spec.ind))
                }
            }
            else {
                for (j in index) {
                     if (rowr[j] == 1 && colr[j] == 1)
                         ret[j] <- ifelse(ind[j] == 1, rbinom(1, 1,
                            sens.ind), 1 - rbinom(1, 1, spec.ind))
                }
            }
        }
        individual <- c(individual, list(indm))
    }
    sq <- rep(1:len, n.col * n.row)
    if (all(colr == 0) && all(rowr == 0))
        return(NULL)
    grd <- data.frame(x = x, col.resp = colr,
        row.resp = rowr, coln = col.groupn, rown = row.groupn,
        arrayn = sq, retest = ret)
    if (ncol(X) > 2)
        for (i in 1:(ncol(X) - 1))
             colnames(grd)[i] <- paste("x", i, sep="")
    list(dframe = grd, ind = individual, prob = as.vector(pijk))
}

