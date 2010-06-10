mp.control <-
function (tol = 0.005, n.gibbs = 1000, n.burnin = 20, 
    maxit = 500, trace = FALSE, time = TRUE) 
{
    if (!is.numeric(tol) || tol <= 0) 
        stop("value of 'tol' must be > 0")
    if (round(n.gibbs) != n.gibbs || n.gibbs <= 0) 
        stop("value of 'n.gibbs' must be a positive integer")
    if (round(n.burnin) != n.burnin || n.burnin <= 0) 
        stop("value of 'n.burnin' must be a positive integer")    
    if (!is.numeric(maxit) || maxit <= 0) 
        stop("maximum number of iterations must be > 0")
    list(tol = tol, n.gibbs = n.gibbs, n.burnin = n.burnin, maxit = maxit, 
        trace = trace, time = time)
}

