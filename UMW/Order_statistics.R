# Order statistics =============================================================

## Order min -------------------------------------------------------------------

or_min<-function(y,alpha,gamma,lambda,n)
{
  min_or<-(alpha*n*y^(-lambda-1)*(-log(y))^gamma*(log(y)*lambda-gamma)*(1-exp(-(alpha*(-log(y))^gamma)/y^lambda))^(n-1)*exp(-(alpha*(-log(y))^gamma)/y^lambda))/log(y)
  return(min_or)
}

## Order max -------------------------------------------------------------------

or_max<-function(y,alpha,gamma,lambda,n)
{
  max_or<-(alpha*n*y^(-lambda-1)*(-log(y))^gamma*(log(y)*lambda-gamma)*exp(-(alpha*n*(-log(y))^gamma)/y^lambda))/log(y)
  return(max_or)
}

## Order r ---------------------------------------------------------------------
or_r<-function(y,alpha,gamma,lambda,r,n)
{
  r_or<-(alpha*y^(-lambda-1)*(-log(y))^gamma*(log(y)*lambda-gamma)*(1-exp(-(alpha*(-log(y))^gamma)/y^lambda))^(n-r)*exp(-(alpha*r*(-log(y))^gamma)/y^lambda))/log(y)
  (factorial(n)/(factorial(r-1)*(factorial(n-r))))*or_max
  return(r_or)
}
