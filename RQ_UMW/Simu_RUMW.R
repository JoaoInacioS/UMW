rm(list=ls())
library(foreach)
source("Gera_RUMW.R")
source("EST_RUMW.R")
source("func_simu_RUMW.R")

#--------- Scenario ---------------------------------------------------------- #

# nc <- c(50,150,300,500)
nc <- c(50,150,300,500) # Required sizes
# tauc <- c(0.1,0.5,0.9)
tauc <- c(0.1,0.5,0.9) # Required Quantile
re<-10100 #Monte Carlo Replicas Simulated
RF<-10000 #Monte Carlo Replicas Required
cen <- 1 # Scenario

{
if(cen==1){theta0<-c(gamma=2.7,lambda=1.8,beta=c(0.2,-0.4,0.5))}
if(cen==2){theta0<-c(gamma=1.5,lambda=2.3,beta=c(0.5,-0.6,0.2))}
}
# link: (logit)
g_lig <- function(mu) {log(mu / (1 - mu))}
ginv_lig <- function(eta) {1 / (1 + exp(-eta))}

# ======= RQ-UMW monte carlo simulation ====================================== #

# start Parallelism ---------------------------------------------------------- #

cores <- parallel::detectCores()
cl <- parallel::makeCluster(cores-1) #not to overload your computer
doSNOW::registerDoSNOW(cl)

# simulation ----------------------------------------------------------------- #

save <- F
nsamplerq <- sim_rumw <- c()
i=1
for (tau in tauc) {
  for (n in nc){
    cat("
tau =",tau,"& n =",n,"----
     ")
    # cov (intercept):
    k<-length(3:length(theta0))-1
    X<-cbind(rep(1,n),matrix(runif(n*k),nrow=n,ncol=k)) # U(0,1)
    #
    nsamplerq[[i]] <- samplerq <- samples_RUMW(theta0 = theta0, X = X, n = n, re = re, tau = tau)
    sim_rumw[[i]] <- suppressWarnings(simulate_RUMW(sample = samplerq, X = X,theta0 = theta0, n = n, 
                                                  tau = tau ,re = re, RF = RF, save = save,intercept = T,
                                                  g_lig = g_lig, ginv_lig = ginv_lig))
    names(nsamplerq)[i] <- names(sim_rumw)[i] <- paste0("n", n)
    i=i+1
  }
cat("-------------------------------------
")
}
# stop paralelo -------------------------------------------------------------- #
foreach::registerDoSEQ()
parallel::stopCluster(cl)



