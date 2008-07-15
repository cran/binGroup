"print.bgtCI" <-
function(x, ...)
{
args<-list(...)
if(is.null(args$digits))
 {digits<-4}
else{digits<-args$digits}
cat("\n The",x$conf.level*100," per cent ", x$method, "-confidence interval is:","\n")
cat("  [",signif(x$conf.int, digits ),"] \n")
cat(" point estimator = ", signif(x$estimate, digits), "\n")
invisible(x)
}

