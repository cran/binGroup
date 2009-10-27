
predict.gt<-function(object, newdata, type=c("link", "response"), se.fit=FALSE, 
                         conf.level=NULL, na.action=na.pass, ...) {

     tt <- terms(object)
     Terms <- delete.response(tt) 
     if (missing(newdata) || is.null(newdata)) {
        m <- model.frame(object)
        newd <- model.matrix(Terms, m)
     }
     else {        
        m <- model.frame(Terms, newdata, na.action = na.action)
        newd <- model.matrix(Terms, m)
     }
     type<-match.arg(type)

     lin.pred<-as.vector(newd%*%object$coefficients)
     link<-object$link
     res<-switch(link, logit=plogis(lin.pred), probit=pnorm(lin.pred), cloglog=1-exp(-exp(lin.pred)))
     if (type=="response") pred<-res 
        else pred<-lin.pred

     if (se.fit) {
        cov<-solve(object$hessian)               
        var.lin.pred<-diag(newd%*%cov%*%t(newd))
        var.res<-switch(link, logit=exp(2*lin.pred)/(1+exp(lin.pred))^4,
             probit=dnorm(lin.pred)^2, cloglog=(exp(-exp(lin.pred))*exp(lin.pred))^2)*var.lin.pred

        if (type=="response") se<-sqrt(var.res) 
           else se<-sqrt(var.lin.pred)

        if (!is.null(conf.level)) {
           alpha<-1-conf.level
           lower<-lin.pred-qnorm(1-alpha/2)*sqrt(var.lin.pred)
           upper<-lin.pred+qnorm(1-alpha/2)*sqrt(var.lin.pred)
           res.lower<-switch(link, logit=plogis(lower), probit=pnorm(lower), cloglog=1-exp(-exp(lower)))
           res.upper<-switch(link, logit=plogis(upper), probit=pnorm(upper), cloglog=1-exp(-exp(upper)))        
           if (type=="response") {
               lwr<-res.lower
               upr<-res.upper
           }
           else {
                lwr<-lower
                upr<-upper
           }
        } 
     }
     if (!is.null(conf.level)) {
        list(fit=pred, se.fit=se, lower = lwr, upper = upr)
     }
     else if (se.fit)
          list(fit=pred, se.fit=se)
     else pred

}
