EM.mp <-
function (col.resp, row.resp, X, coln, rown, sqn, ret, sens,
    spec, linkf, sens.ind, spec.ind, start = NULL, control = mp.control())
{
    if (control$time)
        start.time <- proc.time()
    if (is.null(sens.ind))
        sens.ind <- sens
    if (is.null(spec.ind))
        spec.ind <- spec
    len <- max(sqn)
    diff <- 1
    counts <- 1
    sam <- length(sqn)
    col.groupn <- coln[sqn == 1]
    if (len > 1) {
        for (i in 2:len) {
            temp <- max(col.groupn) + coln[sqn == i]
            col.groupn <- c(col.groupn, temp)
        }
    }
    if (is.null(start)) {
        mod.fit <- try(gtreg.fit(col.resp, X, col.groupn,
            sens, spec, linkf))
        if (class(mod.fit) == "try-error") {
            row.groupn <- rown[sqn == 1]
            if (len > 1) {
                for (i in 2:len) {
                  temp <- max(row.groupn) + rown[sqn == i]
                  row.groupn <- c(row.groupn, temp)
                }
            }
            mod.fit <- gtreg.fit(row.resp, X, row.groupn,
                sens, spec, linkf)
        }
        beta.old <- mod.fit$coefficients
    }
    else beta.old <- start
    extra.loop <- FALSE
    while ((diff > control$tol || extra.loop) && counts <= control$maxit) {
        xib <- X %*% beta.old
        pijk.all <- switch(linkf, logit = plogis(xib),
            probit = pnorm(xib), cloglog = 1 -
                exp(-exp(xib)))
        expect.all <- numeric(0)
        mat2 <- index <- 0        
        erf <- 2 * pijk.all - 1
        for (arrayn in 1:len) {
            index.r <- index.c <- vector("logical", length = sam)
            for (i in 1:sam) {
                if (rown[i] == 1 && sqn[i] == arrayn)
                    index.c[i] <- TRUE
                else index.c[i] <- FALSE
                if (coln[i] == 1 && sqn[i] == arrayn)
                    index.r[i] <- TRUE
                else index.r[i] <- FALSE
            }
            n.row <- max(rown[index.r])
            n.col <- max(coln[index.c])
            rowresp <- row.resp[index.r]
            colresp <- col.resp[index.c]
            index <- max(index) + 1:(n.row * n.col)
            if (!is.null(ret)) {
                re.ind <- na.omit(cbind(coln[sqn == arrayn],
                  rown[sqn == arrayn], ret[sqn == arrayn]))
                re <- ifelse(re.ind[, 3] == 1, sens.ind, 1 -
                  sens.ind)
                re1 <- ifelse(re.ind[, 3] == 0, spec.ind, 1 -
                  spec.ind)
            }
            pijk <- matrix(pijk.all[sqn == arrayn], nrow = n.row)
            a <- ifelse(rowresp == 1, sens, 1 - sens)
            b <- ifelse(colresp == 1, sens, 1 - sens)
            a1 <- ifelse(rowresp == 0, spec, 1 - spec)
            b1 <- ifelse(colresp == 0, spec, 1 - spec)
            mat <- array(NA, c(n.row, n.col, control$n.gibbs))
            y <- matrix(0, nrow = n.row, ncol = n.col)
            for (k in 1:(control$n.gibbs + control$n.burnin)) {
                l <- 1
                for (j in 1:n.col) for (i in 1:n.row) {
                  num <- a[i] * b[j] * pijk[i, j]
                  den.r <- ifelse(sum(y[i, ]) - y[i, j] > 0,
                    a[i], a1[i])
                  den.c <- ifelse(sum(y[, j]) - y[i, j] > 0,
                    b[j], b1[j])
                  den2 <- den.r * den.c * (1 - pijk[i, j])
                  if (!is.null(ret)) {
                    if (l <= length(re) && j == re.ind[l, 1] &&
                      i == re.ind[l, 2]) {
                      num <- num * re[l]
                      den2 <- den2 * re1[l]
                      l <- l + 1
                    }
                  }
                  den <- num + den2
                  if (den != 0) {
                    cond.p <- num/den
                    y[i, j] <- rbinom(1, 1, cond.p)
                  }
                  else y[i, j] <- 0
                }
                if (k > control$n.burnin) {
                  mat[, , k - control$n.burnin] <- y
                  vec <- as.vector(y)
                  for (i1 in index[vec == 1]) for (j1 in index[vec ==
                    1]) {
                    bq <- switch(linkf, logit = 1, probit = 8 *
                      exp(-(xib[i1]^2 + xib[j1]^2)/2)/((1 - erf[i1]^2) *
                      (1 - erf[j1]^2) * pi), cloglog = exp(xib[i1] +
                      xib[j1])/((exp(-exp(xib[i1])) - 1) * (exp(-exp(xib[j1])) -
                      1))) * X[i1, ] %*% t(X[j1, ])
                    mat2 <- mat2 + bq
                  }
                }
            }
            expect.m <- apply(mat, c(1, 2), mean)
            expect <- as.vector(expect.m)
            expect.all <- c(expect.all, expect)
        }
        if (!extra.loop) {
            suppress <- function(w) 
                if(any(grepl("non-integer #successes in a binomial glm", w))) 
                   invokeRestart("muffleWarning")
            mod.fit <- withCallingHandlers(glm.fit(X, expect.all, 
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
    index <- 0
    first <- mat2/control$n.gibbs
    second <- 0
    for (arrayn in 1:len) {
        n.row <- max(rown[sqn == arrayn])
        n.col <- max(coln[sqn == arrayn])
        index <- max(index) + 1:(n.row * n.col)
        expect <- expect.all[index]
        for (i1 in index[expect != 0]) for (j1 in index[expect !=
            0]) {
            coe <- switch(linkf, logit = 1, probit = 8 * exp(-(xib[i1]^2 +
                xib[j1]^2)/2)/((1 - erf[i1]^2) * (1 - erf[j1]^2) *
                pi), cloglog = exp(xib[i1] + xib[j1])/((exp(-exp(xib[i1])) -
                1) * (exp(-exp(xib[j1])) - 1)))
            tim <- expect.all[i1] * expect.all[j1] * coe * X[i1,
                ] %*% t(X[j1, ])
            second <- second + tim
        }
    }
    p2.array <- first - second
    pt1 <- switch(linkf, logit = -exp(xib)/(1 + exp(xib))^2,
        probit = sqrt(2) * xib * exp(-xib^2/2)/(sqrt(pi) * (1 -
            erf)) - 2 * exp(-xib^2)/(pi * (1 - erf)^2), cloglog = -exp(xib))
    pt2 <- switch(linkf, logit = 0, probit = (8 * exp(-xib^2/2) *
        erf + 2 * xib * sqrt(2 * pi) * erf^2 - 2 * xib * sqrt(2 *
        pi)) * exp(-xib^2/2)/((1 + erf)^2 * pi * (1 - erf)^2),
        cloglog = -(exp(xib - exp(xib)) + exp(2 * xib - exp(xib)) -
            exp(xib))/(exp(-exp(xib)) - 1)^2)
    nm <- pt1 + expect.all * pt2
    sign1 <- as.vector(sign(nm))
    nn <- as.vector(sqrt(abs(nm)))
    x2 <- X * nn
    m <- (t(x2) %*% (sign1 * x2))
    H <- -(m + p2.array)
    if (control$time) {
        end.time <- proc.time()
        save.time <- end.time - start.time
        cat("\n Number of minutes running:", round(save.time[3]/60,2), "\n \n")
    }
    list(coefficients = beta.old, hessian = H, df.residual = sam -
        ncol(X), df.null = sam - 1, Gibbs.sample.size = control$n.gibbs,
        counts = count2)
}

