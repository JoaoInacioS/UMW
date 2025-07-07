# Other Functions ==============================================================

## Parallelism -----------------------------------------------------------------

progresso<-function(iterations,sec="1"){
  pb <- progress::progress_bar$new(
    format = paste0("[",sec,"/2]",":percent [:bar] :elapsed | Eta: :eta"),
    total = iterations,    # 100 
    width = 60)
  #foreach: 
  progress <- function(){pb$tick()} 
  opts <- list(progress = progress)
  return(opts)
}

## Function to save simulation -------------------------------------------------

save_cen_RQ<-function(output,RF,tau,n,write){
  l1<-length(output$V1)
  if(l1>=RF){
    total1<-output[1:RF,]
    save(total1,file=paste0("simun_",write,"_",RF,"_",n,"_",tau,".RData"))
  }
  if(l1<RF){
    total1<-output[1:l1,]
    save(total1,file=paste0("simun_",write,"_",l1,"_",n,"_",tau,".RData"))
  }
}

## remove lines with NA --------------------------------------------------------

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

## RQ-UMW function sample replicates -------------------------------------------

samples_RUMW<-function(theta0, X, n, re, tau)
{
  amostra<-bigstatsr::FBM(nrow=n,ncol=re)
  opts<-progresso(iterations = re,sec="1")
  try(foreach(j=1:re, .packages=c("foreach"),
            .combine=rbind,.options.snow = opts,.export = ls(globalenv())) %dopar%{ 
              source("Gera_RUMW.R")
              amostra[,j]<-as.numeric(sample_rumw(n=n,theta0=theta0,X=X,tau=tau))  
      },T)
  amostra<-as.data.frame(amostra[])
  return(amostra)
}

## RQ-UMW monte carlo simulation function --------------------------------------

simulate_RUMW<-function(sample,theta0, X, n, tau, re, RF, intercept=T, save = FALSE,
                        g_lig,ginv_lig)
{
  opts<-progresso(iterations = re,sec="2")
  len_k<-length(theta0)
  k1<-(len_k)*2
  estimate<-bigstatsr::FBM(nrow=re,ncol=k1)
  try(foreach(j=1:re, .packages=c("foreach"),.combine=rbind,.options.snow = opts,
              .export = ls(globalenv())) %dopar% { 
    source("EST_RUMW.R")
    output<-as.data.frame(estimate[])
    output[output==0]<-NA
    l_out<-nrow(which.NA(output))
    if(l_out<=RF-1){
      estimate[j,]<-suppressWarnings(EST_RUMW(y=sample[,j],X=X,tau=tau,intercept = intercept,g_lig=g_lig,ginv_lig=ginv_lig))
    }
    if(l_out>=RF){
    }
  },T)
  output1<-as.data.frame(estimate[])
  output1[output1==0]<-NA
  output1<-which.NA(output1)
  tabela <- data.frame(
    theta0 = theta0,
    Mean = round(colMeans(output1)[1:len_k], 4),
    Bias = round(colMeans(output1)[1:len_k] - theta0, 4),
    RB = round(((colMeans(output1)[1:len_k] - theta0) / theta0) * 100, 4)
  )
  print(tabela)
  # Save estimate
  if(save==T){
    save_cen_RQ(output = output1, RF = RF,  n = n, tau = tau, write = paste0(cen,"_RQUMW"))
  }
  return(list(output=output1,table=tabela))
}





