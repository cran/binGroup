"print.bgtCI" <-
function(x,...)
{
cat("\n The",x$conf.level*100," per cent ", x$method, "-confidence interval is:","\n")
cat("  [",x$conf.int,"] \n")
cat(" point estimator = ", x$estimate, "\n")
invisible(x)
}

