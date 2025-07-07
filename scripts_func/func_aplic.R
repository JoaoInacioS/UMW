# Cleaning function (0.1) ------------------------------------------------------

clean_vector <- function(vector) {
  vector[vector <= 0] <- NA
  vector[vector >= 1] <- NA
  vector <- na.omit(vector)
  return(as.numeric(vector))
}

# Metrics function -------------------------------------------------------------

metrics<-function(y,loglik,par,F_dist,k){
  KS <- function(theta,y,n,fda){
    y <- sort(y)
    aux1 <- c()
    for(i in 1:n){
      aux1[i] <- abs(fda(par=theta,y=y[i])-i/n)
    }   
    return(max(aux1))
  }
  AD <- function(theta,y,n,fda){
    x1 = sort(y)
    u = fda(par=theta,y=x1)
    aux=rep(0,n)
    for(i in 1:n){
      aux[i]=(2*i-1)*log(u[i])+(2*n+1-2*i)*log(1-u[i])
    }
    A2 = -n -(1/n)*sum(aux)
    return(A2)
  }
  cramer <- function(theta,y,n,fda){
    x1 = sort(y)
    u = fda(par=theta,y=x1)
    aux=rep(0,n)
    for(i in 1:n){
      aux[i]=(((2*i-1)/(2*n))-u[i])^2
    }
    W2=sum(aux)+1/(12*n)
    return(W2)
  }
  n<-length(y)
  AIC1 <- -2*loglik+2*k
  BIC1 <- -2*loglik+log(n)*k
  AICc1 <- AIC1 + ((2*k^2)+(2*k))/(n-k-1)
  KS1 <- KS(theta = par,y = y,n = n,fda = F_dist)
  AD1 <- AD(theta = par,y = y,n = n,fda = F_dist)
  CvM1 <- cramer(theta = par,y = y,n = n,fda = F_dist)
  return(list(AIC = AIC1, BIC = BIC1, AICc = AICc1, KS = KS1, AD = AD1, CvM = CvM1))
}

# Summary ----------------------------------------------------------------------

Coef_estim<-function(par,hessian=NULL,f_dist=NULL,y=NULL){
  if(is.null(hessian)){
    if(is.null(f_dist) & is.null(y)){stop("f_dist and y not equal to NULL")}else{
      hessian_dist<-function(par,y,f_dist){
        loglik_func<-function(par,y){sum(log(f_dist(par=par,x = y)))}
        numDeriv::hessian(func=loglik_func,x = par,y=y,method="Richardson")
      }
      hessian<- hessian_dist(par = par,y = y,f_dist = f_dist)
    }}
  stderror<-sqrt(diag(solve(-hessian)))
  #-----------------------z-Values and P-values---------------------------------
  thetaH0<-as.vector(c(rep(0,length(par))))                       
  zstat<-c(abs((par-thetaH0)/stderror))
  pvalues<-2*(1-pnorm(zstat))
  #-----------------------Result Table------------------------------------------
  options(scipen = 1)
  model<-as.data.frame(cbind(round(par,3),round(stderror,3),round(zstat,3),round(pvalues,3)))
  model<-model |> dplyr::mutate(
    nivsig=dplyr::case_when(
      V4 <= 0.001 ~ "***",
      V4 > 0.001 & V4 <= 0.01 ~ "**",
      V4 > 0.01 & V4 <= 0.05 ~ "*",
      V4 > 0.05 & V4 <= 0.1 ~ "·",
      TRUE ~ ""),
    V4=dplyr::case_when(
      V4 < 0.001 ~ paste0("<",0.001),
      TRUE ~ paste0(V4)),
    V2=dplyr::case_when(
      V2 < 0.001 ~ paste0("<",0.001),
      TRUE ~ paste0(V2)))
  colnames(model)<-c("Estimate","Std. Error","z value","Pr(>|z|)","")
  return(model)
}

# Summary RQ-UMW ---------------------------------------------------------------

summary_rqumw<-function(out)
{
  cat(c("============================================================="),fill = TRUE)
  cat(c("Link:",out$link," ","Quantile:",out$quantile," ","Optimization:",out$optimization),fill = TRUE)
  cat(c("============================================================="),fill = TRUE)
  print(out$coef)
  cat(c("---"),fill = TRUE)
  cat(c("Significance code: 0 *** 0.001 ** 0.01 * 0.05 · 0.1 1"),fill = TRUE)
  cat(c("-------------------------------------------------------------"),fill = TRUE)
  cat(c("Loglike:",round(out$loglik,3)," ","AD(p) Resid.:",round(out$metrics$tAD$p.value,3)," ","R2G:",round(out$R2G,3)),fill=TRUE)
  cat(c("AIC:",round(out$metrics$aic,3)," ","BIC:",round(out$metrics$bic,3)," ","AICc:",round(out$metrics$aicc,3)),fill=TRUE)
  cat(c("MSE:",round(out$measures[1],3)," ","RMSE:",round(out$measures[2],3)," ","MAE:",round(out$measures[3],3))," ","MAPE:",round(out$measures[4],3),fill=TRUE)
  cat(c("---"),fill = TRUE)
  cat(c("Number of Iterations in the Optimization Algorithm:",out$outer.iter),fill=TRUE)
  cat((c("Residual standard error: ",round(out$st.res,3)," on ",out$df.residual," degrees of freedom")),fill=TRUE)
  cat("Residuals:",fill=TRUE)
  print(round(summary(as.vector(out$residuals)),3))
  cat(c("============================================================="),fill = TRUE)
}


# Selection of Covariates ======================================================

select_cov<-function(y,X,tau=0.5,ncovn=c(1),writing="",link="logit",
                     criterion=c(AIC,BIC,AICc),intercept=T,save=F){
  names_cov<-c("Cov1","Cov2","Cov3","Cov4","Cov5","Cov6","Cov7","Cov8","Cov9","Cov10")
  models<-c()
  k=1
  for (ncov in ncovn) {
    cat("Number of covariates =",ncov,"(best):
    
")
    models1<-auto_rqUMW(y = y,V = X,tau = tau,ncov = ncov,intercept = intercept,link=link)
    models1<-as.data.frame(models1[])
    models1<-which.NA(models1)
    colnames(models1)<-c(names_cov[1:ncov],"AIC","BIC","AICc","Log_Like","R2G","AD")
    models[[k]]<-models1 |> dplyr::slice_min({{criterion}}, n = 10)
    print(models[[k]][1,])
    k=k+1
    cat("---","
")
  }
  names(models) <- paste0("ncov_",ncovn)
  if(save==T){save(models,file = paste0(writing,"_models_",tau,".RData"))}
  return(models)
}

## auto selection of covariates ------------------------------------------------

auto_rqUMW<-function(y,V,tau,ncov,link="logit",intercept=T)
{
  k<-t(combn(length(V[1,]),ncov))
  nmod<-length(k[,1])
  models<-bigstatsr::as_FBM(matrix(nrow=nmod,ncol=ncov+6))
  opts<-progresso(iterations = nmod)
  foreach(j=1:nmod, .packages=c("foreach"),.combine=rbind
          ,.options.snow = opts,.export = ls(globalenv())) %dopar%{ 
            source("RQ_UMW/EST_RUMW.R")
            source("UMW/EST_semi_closed.R")
            b<-as.numeric(k[j,])
            b1<-na.omit(cbind(y,V[,b]))
            mod<-suppressWarnings(fit_RQUMW(y = b1[,1],X = b1[,2:(ncol(b1))],tau = tau,
                               printmodel=F,link = link,intercept = intercept)) 
            if(class(mod)=="RQ-UMW"){
              if(sum(mod$coef[(1):(ncov+3),4]<0.1,na.rm=TRUE)==ncov+3){
                if(mod$metrics$tAD$p.value>=0.05){
                  models[j,]<-c(b,mod$metrics$aic,mod$metrics$bic,mod$metrics$aicc,
                                mod$loglik,mod$R2G,mod$metrics$tAD$p.value)
                }
              }
            }
          }
  return(models)
}

# Error measures ---------------------------------------------------------------

measures <- function(actual, predicted) {
  mse <- suppressWarnings(mean((actual - predicted)^2))
  rmse <- suppressWarnings(sqrt(mean((actual - predicted)^2)))
  mae <- suppressWarnings(mean(abs(actual - predicted)))
  mape <- suppressWarnings(mean(abs((actual - predicted) / actual)) * 100)
  results <- c(MSE = mse,RMSE = rmse,MAE = mae,MAPE = mape)
  return(results)
}

# R2 generalized coefficient of determination (nagelkerke) ---------------------

R2G<-function(loglik,loglik0,n){
    R2<- 1-exp((-2/n)*((as.numeric((loglik)))-(as.numeric(loglik0))))
    return(c(R2G=R2))
}

# remove lines with NA ---------------------------------------------------------

which.NA<-function(x)
{
  x<-data.frame(x)
  lines<-GLDEX::which.na(x[,1])
  tmp<-length(lines)
  if(tmp>0){
    y<-x[-lines,]
  }
  else{y<-x}
  return(y)
}

# Parallelism ------------------------------------------------------------------

progresso<-function(iterations)
{
  pb <- progress::progress_bar$new(
    format = paste0(":percent [:bar] :elapsed | Eta: :eta"),
    total = iterations,    # 100 
    width = 60)
  #foreach: 
  progress <- function(){pb$tick()} 
  opts <- list(progress = progress)
  return(opts)
}


