#==============================================================================#
# Compiles results, calculating mean, bias, mse, skewness, kurtosis and RB
#==============================================================================#

SIM.COMP<-function(coef,par,n,names.pars=NULL,conf=0.95,arquivo=NULL)
{
  library(knitr)
  library(xtable)
  #----------------------------------------------------------------------------#
  theta0<-par
  thetaH0<-par# as.vector(rep(0,length(theta0)))
  re<-length(coef[,1])
  ne<-length(coef[1,])/2
  coeff<-coef[,1:ne]
  var.coef<-coef[,(ne+1):length(coef[1,])]
  q.norm<-qnorm((1-conf)/2,lower.tail = FALSE)
  #-------Calculation of accuracy measures-------------------------------------#
  vmean<-colMeans(coeff)
  vbias<-vmean-theta0
  vRB<-(((vmean-theta0)/theta0))*100
  vsd<-apply(coeff,2,sd)
  vmse<-vbias^2+vsd^2
  vac<-moments::skewness(coeff)
  vk<-moments::kurtosis(coeff)
  #-----------IC e TH----------------------------------------------------------#
  # contIC<-matrix(0,ncol=length(theta0),nrow=1)
  contIC<-bigstatsr::FBM(nrow=re,ncol=length(theta0))
  opts<-progresso(iterations = re)
  foreach(i=1:re, .packages=c("foreach"),.combine=rbind,.options.snow = opts) %dopar%{ 
            LS<-coeff[i,]+q.norm*sqrt(var.coef[i,])
            LI<-coeff[i,]-q.norm*sqrt(var.coef[i,])
            sum_1<-(LI<=theta0 & theta0<=LS)
            contIC[i,]<-as.numeric(ifelse(sum_1==TRUE,1,0))
  }
  tc<-as.numeric(colSums(as.data.frame(contIC[])))/re
  results<-matrix(0,nrow=8,ncol=length(theta0))
  results<-round(rbind(theta0,vmean,vbias,vRB,vsd,vmse,vac,vk,tc),3)
  if(!is.null(names.pars)){colnames(results)<-names.pars}
  P<-c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","P11")
  if(is.null(names.pars)){colnames(results)<-P[1:ne]}
  NIC<-paste0("CR",conf*100,"%")
  rownames(results)<-c("Theta0","Mean","Bias","RB","SD","MSE","AC","K",NIC)
  if(!is.null(arquivo)){
    #-----------saves the Results Table to the file----------------------------#
    # write(c("#---------------------------"),file=arquivo,append=TRUE,ncolumns=1)
    # print(xtable(results,type="latex",digits = 4),file=arquivo,append=TRUE)
    write(c(paste0("# n = ",n," ----")),file=arquivo,append=TRUE,ncolumns=1)
    print(xtable(t(results),type="latex",digits = 4),file=arquivo,append=TRUE)
  }else{
    return(print(kable(t(results))))
  }
}

# Parallelism -----------------------------------------------------------------#

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


# UMW ==========================================================================

scenarios<-c()
scenarios[[1]]<-c(alpha=0.7,gamma=1.3,lambda=0.5)
scenarios[[2]]<-c(alpha=0.3,gamma=0.8,lambda=1.2)
scenarios[[3]]<-c(alpha=1.3,gamma=1.1,lambda=0.6)
scenarios[[4]]<-c(alpha=0.5,gamma=0.9,lambda=0.8)
ncen<-c(1,2,3,4) # required scenarios
nc<-c(40,80,120,160,200) # required n
replicas<-10000
#
cores <- parallel::detectCores()
cl <- parallel::makeCluster(cores-1) #not to overload your computer
doSNOW::registerDoSNOW(cl)
library(foreach)
#
for (cen in ncen) {
  theta0<-scenarios[[cen]]
  cat("
====================================>",cen,"<=====================================
")
  for (n in nc) {
    cat("# n =",n,":
")
    load(paste0("./UMW/simun_",cen,"_UMW_",replicas,"_",n,".RData"))
    (SIM.COMP(coef = total1,n=n,par = theta0,arquivo = NULL,names.pars = c("alpha","gamma","lambda")))
    cat("---
")
  }
}
#
foreach::registerDoSEQ()
parallel::stopCluster(cl)


# RQ-UMW =======================================================================

scenarios<-c()
scenarios[[1]]<-c(gamma=2.7,lambda=1.8,beta=c(0.2,-0.4,0.5))
scenarios[[2]]<-c(gamma=1.5,lambda=2.3,beta=c(0.5,-0.6,0.2))
ncen<-c(1,2) # required scenarios
nc<-c(50,150,300,500) # required n
tauc<-c(0.1,0.5,0.9) # required quantile
replicas<-10000
#
cores <- parallel::detectCores()
cl <- parallel::makeCluster(cores-1) #not to overload your computer
doSNOW::registerDoSNOW(cl)
library(foreach)
#
for (cen in ncen) {
  theta0<-scenarios[[cen]]
  cat("
====================================>",cen,"<=====================================
")
  for (tau in tauc) {
    for (n in nc) {
      cat("# tau =",tau,"& n =",n,":
")
      load(paste0("./RQUMW/simun_",cen,"_RQUMW_",replicas,"_",n,"_",tau,".RData"))
      (SIM.COMP(coef = total1,n=n,par = theta0,arquivo = NULL,names.pars = c("gamma","lambda",paste0("beta",0:2))))
      cat("---
")
    }
    cat("-----------------------------------------------------------------------------
")
  }
}
#
foreach::registerDoSEQ()
parallel::stopCluster(cl)

