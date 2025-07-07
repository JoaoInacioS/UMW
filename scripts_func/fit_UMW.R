# Fit UMW ======================================================================
source("UMW/Est_semi_closed.R")
source("scripts_func/func_aplic.R")
#------------------------------------------------------------------------------#
fit_UMW<-function(y){
  out<-c()
  n_x <- length(y)
  if(any(is.na(y))) {
    num_na <- sum(is.na(y))
    out$num_na <- num_na
  }else {
    out$num_na <- 0
  }
  if(any(y <= 0) | any(y >= 1) | any(is.na(y))){
    y <- clean_vector(vector = y)
    n_x2 <- length(y)
    num_omit <- n_x - n_x2 - out$num_na
    out$num_omit <- num_omit
  }else {
    out$num_omit <- 0
  }
  F_umw<-function(par,y){
    alpha<-par[1]
    gamma<-par[2]
    lambda<-par[3]
    exp(-(alpha*(-log(y))^gamma)/y^lambda)
  }
  #
  out$summary <- c(n = length(y),round(c(summary(y)[c(1, 3, 4, 6)],
                   sd = sd(y),AC = moments::skewness(y),K = moments::kurtosis(y),
                   AD = nortest::ad.test(y)$p.value), 3))
  mod1<-suppressWarnings(Est_UMW(y=y,applic = TRUE))
  out$metrics<-metrics(y=y,loglik = mod1$value,par = mod1$par,F_dist = F_umw,k=3)
  out$coef<-Coef_estim(par = mod1$par,hessian = mod1$hessian)
  out$convergence<-mod1$convergence==0
  #
  cat("
Summary:
")
  print(out$summary)
  cat("
")
  cat("
Coefficients:
")
  print(out$coef)

  cat("
Metrics:
")
  tab<-round(as.numeric(out$metrics),3)
  names(tab)<-c("AIC","BIC","AICc","KS","AD","CvM")
  print(tab)
  out2<-out
}


