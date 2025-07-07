# Density and Cumulative =======================================================

## Probability Density Function UMW ----
dUMW<-function(y,alpha,gamma,lambda)
{
  f<- (alpha/(log(y)*y^(lambda+1)))*((-log(y))^gamma)*(lambda*log(y)-gamma)*exp(-alpha*((-log(y))^gamma)*y^(-lambda))
  return(f)
}

## Cumulative Function UMW ----
pUMW<-function(y,alpha,gamma,lambda)
{
  f<- exp(-alpha*((-log(y))^gamma)*(y^(-lambda)))
  return(f)
}


# generate UMW observations ====================================================

umw<-function(n,theta0)
{
  alpha<-theta0[1]
  gamma<-theta0[2]
  lambda<-theta0[3]
  eq<-function(x,parms)
    {((-log(x))^parms[3])*(x^(-parms[4]))+(1/parms[2])*log(parms[1])}
  y<-c()
  j<-1
  while(j<=n){
    u<-runif(1)
    parms<-c(u,alpha,gamma,lambda)
    tmp<-try(suppressWarnings(stats::uniroot(f=eq,lower=0.0001,upper=0.9999,parms=parms)),T)
    if(class(tmp)=="list"){
      y[j]<-tmp$root
      j<-j+1
    }
  }
  return(y)
}

