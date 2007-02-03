"print.bgtTest" <-
function(x,...)

{
if (x$alternative=="two.sided"){alt.hyp="true proportion is not equal to "}
if (x$alternative=="less"){alt.hyp="true proportion is less than "}
if (x$alternative=="greater"){alt.hyp="true proportion is greater than "}

cat("\n Alternative hypothesis: ",alt.hyp,x$p.hyp,"\n")
cat(" p-value = ", x$p.val,"\n")
cat(" point estimator = ", x$estimate,"\n")
invisible(x)
}

