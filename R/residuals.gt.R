residuals.gt <-
function(object, type = c("deviance", "pearson", 
    "response"), ...) {

    type <- match.arg(type)
    r <- object$residuals
    zhat<-object$fitted.group.values
    z<-object$z
    res<-switch(type, response= r, pearson=r/sqrt(zhat*(1-zhat)),
           deviance=sqrt(-2*(z*log(zhat)+(1-z)*log(1-zhat)))*sign(z-zhat))    
    res

}

