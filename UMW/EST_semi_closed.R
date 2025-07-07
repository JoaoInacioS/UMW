
# Est UMW ======================================================================

Est_UMW<-function(y,applic = F)
{
  lb <- c(rep(0.001,2))
  ub <- c(rep(Inf,2))
  start.theta<-c(1,1)
  mod1<-try(optim(par=start.theta,fn=llike_UMW,lower=lb,upper=ub,y=y,
                  method="L-BFGS-B",gr=vscore_UMW,hessian=F,m.optim=1.0,
                  control=list(fnscale=-1,maxit=2000)),T)
  alpha<-length(y)/sum((-log(y))^mod1$par[1]/y^mod1$par[2])
  mod1$par<-c(alpha,mod1$par)
  mod1$hessian<-hessian_UMW(theta0=mod1$par,y=y)
  tmp2<-test.fun(mod1)
  if(class(tmp2)=="numeric"){
    if(applic==F){tmp2}else{return(mod1)}
  }else{if(applic==F){rep(NA,k1)}else{return(cat("ALGORITHM DID NOT CONVERGE!"))}}
}

## Log-vero UMW ----------------------------------------------------------------

llike_UMW<-function(theta0,y,m.optim = 1)
{
  gamma<-theta0[1]
  lambda<-theta0[2]
  alpha<-length(y)/sum((-log(y))^gamma/y^lambda)
  log_like<-sum(log(alpha)-(lambda+1)*log(y)+(gamma-1)*log(-log(y))+log(gamma-lambda*log(y))-
                    alpha*((-log(y))^gamma)*y^(-lambda))
  if(m.optim==-1){return(-log_like)} 
  if(m.optim==1){return(log_like)}
}

## Score Function MLE ----------------------------------------------------------

vscore_UMW<-function(theta0,y,m.optim = 1)
{
  gamma<-theta0[1]
  lambda<-theta0[2]
  alpha<-length(y)/sum((-log(y))^gamma/y^lambda)
  # alpha
  wt<- (1/alpha)-((-log(y))^gamma/y^lambda)
  # gamma
  rt<- 1/(gamma-log(y)*lambda)-(alpha*(-log(y))^gamma*log(-log(y)))/y^lambda+log(-log(y))
  # lambda 
  st <- -log(y)/(gamma-log(y)*lambda)+(alpha*(-log(y))^gamma*log(y))/y^lambda-log(y)
  Ugamma<-sum(rt); Ulambda<-sum(st)
  vetor_score<-c(Ugamma,Ulambda)
  return(vetor_score)
}

## Hessian UMW -----------------------------------------------------------------

hessian_UMW<-function(theta0,y)
{
  alpha<-theta0[1]
  gamma<-theta0[2]
  lambda<-theta0[3]
  #
  n=length(y)
  I_n<-(rep(1,n))
  d_aa<-(-n/alpha^2)
  d_gg<-(-1/(gamma-log(y)*lambda)^2-(alpha*(-log(y))^gamma*log(-log(y))^2)/y^lambda)
  d_ll<-(-log(y)^2/(gamma-log(y)*lambda)^2-(alpha*(-log(y))^gamma*log(y)^2)/y^lambda)
  d_ag<-(-((-log(y))^gamma*log(-log(y)))/y^lambda)
  d_al<-(((-log(y))^gamma*log(y))/y^lambda)
  d_gl<-(log(y)/(gamma-log(y)*lambda)^2+(alpha*(-log(y))^gamma*log(y)*log(-log(y)))/y^lambda)
  #
  J_aa<-d_aa
  J_gg<-t(d_gg)%*%I_n
  J_ll<-t(d_ll)%*%I_n
  J_ag<-t(d_ag)%*%I_n
  J_al<-t(d_al)%*%I_n
  J_gl<-t(d_gl)%*%I_n
  #
  (matrix(c(J_aa,J_ag,J_al,J_ag,J_gg,J_gl,J_al,J_gl,J_ll),nrow = 3,byrow = T))
}

# Other Functions ==============================================================

## all positions of the vector are positive ------------------------------------

is.positive<-function(a)
{
  k<-length(a)
  tmp<-sum(a>0)
  return(k==tmp)
}

## Convergence test ------------------------------------------------------------

test.fun<-function(object)
{
  if(class(object)=="list"){
    if(object$convergence==0){
      parameters<-try(object$par,T)
      hess<-try(object$hessian,T)
      var.coef<-try(diag(solve(-hess)),T)
      if(is.numeric(parameters)==TRUE){
        if(is.numeric(var.coef)==TRUE){
          if(is.positive(var.coef)==TRUE){
            z<-c(parameters,var.coef)
            return(z)
          }else{return(FALSE)}
        }else{return(FALSE)}
      }else{return(FALSE)}
    }else{return(FALSE)}
  }else{return(FALSE)}
}








