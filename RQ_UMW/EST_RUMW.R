# Est RQ-UMW ===================================================================

EST_RUMW<-function(y,X,tau=0.5,applic=F,intercept=T,
                   g_lig,ginv_lig)
{
  #----------Start Point an Limits For the Parameters--------------------------#
  ynew<-g_lig(y)
  # lm: Ordinary Least Squares
  Y<-(X[,-1])
  df <- data.frame(ynew, Y)
  if(intercept==T){mod.ols<-try(lm(ynew~.,data=df),T)
  }else{mod.ols<-try(lm(ynew~.-1,data=df),T)}
  startbeta<-mod.ols$coefficients
  #
  lb <- c(rep(0.001,2),rep(-Inf,length(startbeta)))
  ub <- c(rep(Inf,length(startbeta)+2))
  start.theta<-c(1,1,startbeta)
  mod<-try(optim(par=start.theta,fn=llike_RUMW,gr=vscore_RUMW,lower=lb,upper=ub,y=y,X=X,
                 method="L-BFGS-B",tau=tau,hessian=F,m.optim=1.0,ginv_lig=ginv_lig,g_lig=g_lig,
                 control=list(fnscale=-1,maxit=2000)),T)
  mod$hessian<-hessian_RUMW(theta0=mod$par,y=y,X=X,tau=tau,ginv_lig=ginv_lig,g_lig=g_lig)
  tmp2<-test.fun(mod)
  if(class(tmp2)=="numeric"){
    if(applic==F){tmp2}else{return(mod)}
  }else{if(applic==F){rep(NA,k1)}else{return(cat("ALGORITHM DID NOT CONVERGE!"))}}
}

## Log-vero RQ-UMW -------------------------------------------------------------

llike_RUMW<-function(theta0,y,X,tau,m.optim,ginv_lig,g_lig)
{
  beta<-theta0[3:length(theta0)]
  mut<-ginv_lig(as.vector(X%*%beta))
  gamma<-theta0[1]
  lambda<-theta0[2]
  alpha <- -((mut^lambda)*log(tau))/((-log(mut)))^(gamma)  # reparameterization
  log_like<-sum(log(alpha)-(lambda+1)*log(y)+(gamma-1)*log(-log(y))+log(gamma-lambda*log(y))-
                alpha*((-log(y))^gamma)*y^(-lambda))
  if(m.optim==-1){return(-log_like)} 
  if(m.optim==1){return(log_like)}
}

## Scoring function RQ-UMW -----------------------------------------------------

vscore_RUMW<-function(theta0,y,X,tau,g_lig,ginv_lig,m.optim = 1,vsmatrix=F)
{
  D1_num <- Deriv::Deriv(g_lig)
  beta<-theta0[3:length(theta0)]
  mut<-ginv_lig(as.vector(X%*%beta))
  gamma<-theta0[1]
  lambda<-theta0[2]
  # gamma
  rt<- 1/(gamma-log(y)*lambda)+log(-log(y))-log(-log(mut))+
  (mut^lambda*log(tau)*(-log(y))^gamma*(log(-log(y))-log(-log(mut))))/((-log(mut))^gamma*y^lambda)
  # lambda 
  st <- -log(y)/(gamma-log(y)*lambda)-log(y)+log(mut)+
  (mut^lambda*log(tau)*(log(mut)-log(y))*(-log(y))^gamma)/((-log(mut))^gamma*y^lambda)
  # mu
  wt <- ((y^lambda*(-log(mut))^gamma+log(tau)*(-log(y))^gamma*mut^lambda)*(lambda*log(mut)-gamma))/(y^lambda*mut*(-log(mut))^gamma*log(mut))
  if(vsmatrix == F){
    Ugamma<-sum(rt); Ulambda<-sum(st); Ubeta <- as.vector(t(X) %*% (diag(1/D1_num(mut))) %*% wt)
    vetor_score<-c(Ugamma,Ulambda,Ubeta)
    return(vetor_score)
  }else{
    Ubeta_m <- X * (mut*(1-mut) * wt)
    m_score<-as.matrix(cbind(rt,st,Ubeta_m))
    return(m_score)
  }
}

## Hessian RQ-UMW --------------------------------------------------------------

hessian_RUMW<-function(theta0,y,X,g_lig,ginv_lig,tau,m.optim = 1)
{
  D1_num <- Deriv::Deriv(g_lig)
  D2_num <- Deriv::Deriv(D1_num)
  n<-length(y)
  beta<-theta0[3:length(theta0)]
  mut<-ginv_lig(as.vector(X%*%beta))
  gamma<-theta0[1]
  lambda<-theta0[2]
  T1<-(diag(1/D1_num(mut)))
  T2<-(diag(-(D2_num(mut)*(1/(D1_num(mut))^2))))
  A_1 <- mut^lambda*log(tau)*(-log(y))^gamma
  A_2 <- log(-log(y))-log(-log(mut))
  #
  d_GG<-((A_1*A_2^2)/((-log(mut))^gamma*y^lambda)-1/(gamma-log(y)*lambda)^2)
  d_LL<-((A_1*(log(mut)-log(y))^2/((-log(mut))^gamma*y^lambda))-log(y)^2/(gamma-log(y)*lambda)^2)
  d_GL<-(log(y)/(gamma-log(y)*lambda)^2+(A_1*A_2*(log(mut)-log(y)))/((-log(mut))^gamma*y^lambda))
  d_Gb<-((-((A_1*(A_2*(gamma-log(mut)*lambda)+1))/((-log(mut))^gamma*y^lambda))-1)/(mut*log(mut)))
  d_Lb<-(((-log(mut))^gamma*y^lambda+A_1)/(mut*(-log(mut))^gamma*y^lambda)-(A_1*(log(y)-log(mut))*(log(mut)*lambda-gamma))/(mut*(-log(mut))^gamma*log(mut)*y^lambda))
  d_bb<-(gamma/(mut^2*log(mut)^2)-(log(mut)*gamma*(A_1*(2*lambda-1)-(-log(mut))^gamma*y^lambda)-A_1*gamma*(gamma+1))/(mut^2*(-log(mut))^gamma*log(mut)^2*y^lambda)) -
    ((A_1*(1-lambda)+(-log(mut))^gamma*y^lambda)*lambda)/(mut^2*(-log(mut))^gamma*y^lambda)
  wt <- ((y^lambda*(-log(mut))^gamma+log(tau)*(-log(y))^gamma*mut^lambda)*(lambda*log(mut)-gamma))/(y^lambda*mut*(-log(mut))^gamma*log(mut))
  #
  I_n<-(rep(1,n))
  J_GG<-(d_GG)%*%I_n
  J_LL<-t(d_LL)%*%I_n
  J_GL<-t(d_GL)%*%I_n
  J_Gb<-t(X)%*%T1%*%d_Gb
  J_Lb<-t(X)%*%T1%*%d_Lb
  T3<-diag(c((d_bb)%*%T1+(wt)%*%T2))
  J_bb<-t(X)%*%T3%*%T1%*%X
  #
  P1<-(matrix(c(J_GG,J_GL,t(J_Gb),J_GL,J_LL,t(J_Lb)),nrow = 2,ncol = length(beta)+2,byrow = T))
  P2<-cbind(J_Gb,J_Lb,J_bb)
  hessian<-rbind(P1,P2)
  return(hessian)
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
