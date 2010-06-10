EM.control <-
function (tol = 0.001, maxit = 1000, trace = FALSE) 
{
    if (!is.numeric(tol) || tol <= 0) 
        stop("value of 'tol' must be > 0")
    if (!is.numeric(maxit) || maxit <= 0) 
        stop("maximum number of iterations must be > 0")
    list(tol = tol, maxit = maxit, trace = trace)
}

