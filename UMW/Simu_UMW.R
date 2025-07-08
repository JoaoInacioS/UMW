rm(list=ls())
library(foreach)
source("GERA_UMW.R")
source("EST_semi_closed.R")
source("func_simu_UMW.R")

#--------- Scenario ---------------------------------------------------------- #

# nc <- c(40,80,120,160,200)
nc <- c(40,80,120,160,200) # Required sizes
re<-10100 #Monte Carlo Replicas Simulated
RF<-10000 #Monte Carlo Replicas Required
cen <- 1 # Scenario

{
if(cen==1){theta0<-c(alpha=0.7,gamma=1.3,lambda=0.5)} # Scenario 1
if(cen==2){theta0<-c(alpha=0.3,gamma=0.8,lambda=1.2)} # Scenario 2
if(cen==3){theta0<-c(alpha=1.3,gamma=1.1,lambda=0.6)} # Scenario 3
if(cen==4){theta0<-c(alpha=0.5,gamma=0.9,lambda=0.8)} # Scenario 4
}

# ======= UMW monte carlo simulation ========================================= #

# start Parallelism ---------------------------------------------------------- #

cores <- parallel::detectCores()
cl <- parallel::makeCluster(cores-1) #not to overload your computer
doSNOW::registerDoSNOW(cl)

# simulation ----------------------------------------------------------------- #

save <- F
nsample <- sim_umw <- c()
i=1
for (n in nc){
  cat("
n =",n,"----
   ")
  nsample[[i]] <- sample1 <- sample_UMW(theta0 = theta0, n = n, re = re)
  sim_umw[[i]] <- suppressWarnings(simulate_UMW(sample = sample1,theta0 = theta0, n = n, 
                               re = re, RF = RF, save = save))
  names(nsample)[i] <- names(sim_umw)[i] <- paste0("n", n)
  i=i+1
}

# stop paralelo -------------------------------------------------------------- #
foreach::registerDoSEQ()
parallel::stopCluster(cl)




