# Other Functions ==============================================================

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

## Parallelism -----------------------------------------------------------------

progresso<-function(iterations,sec="1")
{
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

save_cen<-function(output,RF,n,write="")
{
  l1<-length(output$V1)
  if(l1>=RF){
    total1<-output[1:RF,]
    save(total1,file=paste0("simun_",write,"_",RF,"_",n,".RData"))
  }
  if(l1<RF){
    total1<-output[1:l1,]
    save(total1,file=paste0("simun_",write,"_",l1,"_",n,".RData"))
  }
}

## sample replicate function ----------------------------------------------------

sample_UMW<-function(theta0=c(0.7,1.3,0.5), n = 40, re = 100)
{
  amostra<-bigstatsr::FBM(nrow=n,ncol=re)
  opts<-progresso(iterations = re,sec = "1")
  try(foreach(j=1:re, .packages=c("foreach"),.combine=rbind,.options.snow = opts) %dopar%{ 
    source("Gera_UMW.R")
    amostra[,j]<-as.numeric(umw(n=n,theta0=theta0))  
  },T)
  amostra<-as.data.frame(amostra[])
  return(amostra)
}


## UMW monte carlo simulation function -----------------------------------------

simulate_UMW<-function(sample,theta0=c(0.7,1.3,0.5), n = 40, re = 100, RF = 100, save = FALSE)
{
  opts<-progresso(iterations = re,sec = "2")
  k1<-(length(theta0))*2
  estimate<-bigstatsr::FBM(nrow=re,ncol=k1)
  try(foreach(k=1:re, .packages=c("foreach"),.combine=rbind,.options.snow = opts,.export = ls(globalenv())) %dopar% { 
    source("Est_semi_closed.R")
    output<-as.data.frame(estimate[])
    output[output==0]<-NA
    l_out<-nrow(which.NA(output))
    if(l_out<RF){
      estimate[k,]<-suppressWarnings(Est_UMW(y = sample[,k]))
    }
    if(l_out>=RF){
    }
  },T)
  output1<-as.data.frame(estimate[])
  output1[output1==0]<-NA
  output1<-which.NA(output1)
  tabela <- data.frame(
    theta0 = theta0,
    Mean = round(colMeans(output1)[1:3], 4),
    Bias = round(colMeans(output1)[1:3] - theta0, 4),
    RB = round(((colMeans(output1)[1:3] - theta0) / theta0) * 100, 4)
  )
  print(tabela)
  # Save estimate
  if(save==T){
    save_cen(output = output1, RF = RF,  n = n, write = paste0(cen,"_UMW"))
  }
  return(list(output=output1,table=tabela))
}
