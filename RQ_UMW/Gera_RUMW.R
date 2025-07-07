# Generate reparameterized sample RQ-UMW =======================================

sample_rumw<-function(n,theta0,X,tau=0.5)
{
  beta<-theta0[3:length(theta0)]
  mut<-(1/(1+exp(-X%*%beta)))
  gamma<-theta0[1]
  lambda<-theta0[2]
  alpha <- -((mut^lambda)*log(tau))/((-log(mut)))^(gamma)  # reparameterization
  xt<-rumw(n=n,alpha=alpha,gamma=gamma,lambda=lambda)
  return(as.vector(xt))
}

## generate UMW observations (regression) --------------------------------------

rumw<-function(n,alpha,gamma,lambda)
{
  eq<-function(x,parms)
  {((-log(x))^parms[3])*(x^(-parms[4]))+(1/parms[2])*log(parms[1])}
  y<-c()
  j<-1
  while(j<=n){
    u<-runif(1)
    parms<-c(u,alpha[j],gamma,lambda)
    tmp<-try(suppressWarnings(stats::uniroot(f=eq,lower=0.0001,upper=0.9999,parms=parms)),T)
    if(class(tmp)=="list"){
      y[j]<-tmp$root
      j<-j+1
    }
  }
  return(y)
}





