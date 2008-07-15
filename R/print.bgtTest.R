"print.bgtTest" <-
function(x,...)

{
args<-list(...)
if(is.null(args$digits))
 {digits<-4}
else{digits<-args$digits}

if (x$alternative=="two.sided"){alt.hyp="true proportion is not equal to "}
if (x$alternative=="less"){alt.hyp="true proportion is less than "}
if (x$alternative=="greater"){alt.hyp="true proportion is greater than "}

cat("\n Alternative hypothesis: ",alt.hyp,x$p.hyp,"\n")
cat(" p-value = ",signif( x$p.val, digits),"\n")
cat(" point estimator = ", signif( x$estimate, digits),"\n")
invisible(x)
}

