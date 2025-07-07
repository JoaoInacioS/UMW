# Fit UMW ======================================================================
source("RQ_UMW/EST_RUMW.R")
source("UMW/EST_semi_closed.R")
source("scripts_func/func_aplic.R")
#------------------------------------------------------------------------------#
fit_RQUMW<-function(y,X,tau,link="logit",intercept=T,printmodel=T)
{
  #----------------------------------------------------------------------------#
  y<-as.numeric(na.omit(y))
  if(min(y) < 0 & max(y) > 0){stop("OUT OF RANGE!")}
  #-----------------------Output-----------------------------------------------#
  out<-c()
  out$X=X
  out$y<-y
  out$quantile<-tau
  out$n<-length(y)
  out$optimization<-"L-BFGS-B"
  if(class(link)=="function"){out$link<-"unknown"}else{out$link<-link}
  class(out)<-c("RQ-UMW")
  #-----------------------------intercept--------------------------------------#
  if(intercept==T){X<-as.matrix(cbind(`(intercept)`=1,X))}else{X<-as.matrix(X)}
  #-----------------------Choice of Link Function------------------------------#
  finv_num <- function(link_function, eta_value, start_value = rep(0.5,length(eta_value))) {
    inverse_function <- function(mu) {
      link_function(mu) - eta_value
    }
    result <- suppressWarnings(nleqslv::nleqslv(start_value, inverse_function)$x)
    return(result)
  }
  if(class(link)=="character"){
    if(link == "logit"){
    g_lig <- function(mu) {log(mu / (1 - mu))}
    ginv_lig <- function(eta) {1 / (1 + exp(-eta))}
    }
    if(link=="probit"){
    g_lig <- function(mu) {qnorm(mu)}
    ginv_lig <- function(eta) {pnorm(eta)}
    }
    if(link=="cauchit"){
    g_lig <- function(mu) {tan(pi * (mu - 0.5))}
    ginv_lig <- function(eta) {0.5 + (atan(eta) / pi)}
    }
    if(link=="cloglog"){
    g_lig <- function(mu) {log(-log(1 - mu))}
    ginv_lig <- function(eta) { 1 - exp(-exp(eta))}
    }
    if(link=="loglog"){
    g_lig <- function(mu) {log(-log(mu))}
    ginv_lig <- function(eta) {exp(-exp(eta))}
    }
  }
  if(class(link)=="function"){
    g_lig <- link
    ginv_lig <- function(eta) {finv_num(g_lig,eta_value = eta)}
  }
  out$g_lig<-g_lig
  out$ginv_lig<-ginv_lig
  if(class(link)!="function" & class(link)!="character"){
    stop("Return the link function as a 'function' class or using the definitions 
         'logit', 'probit', 'cauchit', 'loglog', 'cloglog'.")
  }
  D1_num <- Deriv::Deriv(g_lig)
  D2_num <- Deriv::Deriv(D1_num)
  #-----------------------colnames(X)------------------------------------------#
  beta0<-c(expression(beta[0]),expression(beta[1]),expression(beta[2]),expression(beta[3]),
           expression(beta[4]),expression(beta[5]),expression(beta[6]),expression(beta[7]))  
  tmpnames<-c(beta0[1:ncol(as.matrix(X))])
  if(!is.null(colnames(X))){tmpnames<-colnames(X)}
  #-------------------------Estimate-------------------------------------------#
  tmp1<-suppressWarnings(EST_RUMW(y=y,X=X,tau=tau,applic = TRUE,intercept = intercept,
                                  g_lig = g_lig, ginv_lig = ginv_lig))
  out$loglik<-tmp1$value
  out$convergence<-tmp1$convergence==0
  pars<-tmp1$par
  out$outer.iter<-as.numeric(tmp1$counts[1])
  out$hessian<- tmp1$hessian
  #-----------------------Initial Part of Function Output----------------------#
  beta<-pars[3:length(pars)]
  zetahat<-X%*%as.matrix(beta)
  muhat<-ginv_lig(zetahat)
  gamma<-pars[1]
  lambda<-pars[2]
  names(pars)<-c("gamma","lambda",tmpnames)
  alpha <- -((muhat^lambda)*log(tau))/((-log(muhat)))^(gamma)  # reparameterization
  out$alpha<-alpha
  out$fitted<-as.vector(muhat)
  out$zetahat<-zetahat
  out$qrnames<-names(pars)
  out$pars<-pars
  #-------------------------Estimate0------------------------------------------#
  tmp0<-suppressWarnings(Est_UMW(y=y,applic = T))
  out$loglik0<-tmp0$value
  #-----------------------R^2--------------------------------------------------- 
  out$R2G<-1-exp((-2/out$n)*((out$loglik)-(out$loglik0)))
  #-----------------Cumulative Distribution Function---------------------------#
  peged<-function(x,alpha,gamma,lambda)
  {
    fpa<-(exp(-(alpha*(-log(y))^gamma)/y^lambda))
    return(fpa)
  }
  #------Quantile Residue------------------------------------------------------#
  quantresid<-qnorm(peged(x = y, alpha = alpha, gamma = gamma,lambda = lambda))
  out$residuals<-as.numeric(quantresid)
  out$df.residual<-out$n-length(coef)
  out$st.res<-sum((as.vector(out$residuals)-mean(as.vector(out$residuals)))^2)/out$df.residual
  out$coef<-Coef_estim(par = out$pars,hessian = out$hessian)
  #-----------------------Model Selection Criteria-----------------------------#
  out$metrics$aic <- -2*out$loglik+2*(2+length(beta))
  out$metrics$bic <- -2*out$loglik+log(out$n)*(2+length(beta))
  out$metrics$aicc <- out$metrics$aic + (((2*(2+length(beta))^2)+(2*(2+length(beta))))/(out$n-(2+length(beta))-1))
  out$metrics$tAD <- nortest::ad.test(out$residuals)
  #----------------------------Error measures----------------------------------#
  out$measures <- measures(actual = out$y,predicted = out$fitted)
  #-----------------------Result Table-----------------------------------------#
  if(printmodel==TRUE){
    summary_rqumw(out)
  }  
  out2<-out
}
