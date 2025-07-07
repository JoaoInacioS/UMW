rm(list=ls())
source("scripts_func/fit_UMW.R")
source("scripts_func/fit_RQUMW.R")
source("scripts_func/func_aplic.R")

# Empirical Applications =======================================================

## SDG 17.3 --------------------------------------------------------------------

# total municipal revenues collected of the municipalities of 
# Rio Grande do Sul (RS) state in the year 2021

data1 <- read.table("data_set/SDG17_3_P_RC_TRB_2021.txt",header=T)

y1<-data1$SDG17_3_P_RC_TRB_2021
Aj<-fit_UMW(y = y1)


## Useful Volume ---------------------------------------------------------------

Useful<-read.table("data_set/useful.txt", header = TRUE, sep = "\t")

### South ---- 

AjS1<-fit_UMW(Useful$BARRA.GRANDE)
AjS2<-fit_UMW(Useful$CAMPOS.NOVOS) 
AjS3<-fit_UMW(Useful$G..B..MUNHOZ)
AjS4<-fit_UMW(Useful$G..P..SOUZA)
AjS5<-fit_UMW(Useful$MACHADINHO)
AjS6<-fit_UMW(Useful$MAUA)
AjS7<-fit_UMW(Useful$PASSO.FUNDO)
AjS8<-fit_UMW(Useful$SALTO.SANTIAGO)
AjS9<-fit_UMW(Useful$SANTA.CLARA.PR)
AjS10<-fit_UMW(Useful$SEGREDO)

### Southeat/Midwest ---- 

AjSeMw1<-fit_UMW(Useful$AVERMELHA)
AjSeMw2<-fit_UMW(Useful$B..BONITA)
AjSeMw3<-fit_UMW(Useful$CAPIVARA)
AjSeMw4<-fit_UMW(Useful$CHAVANTES)
AjSeMw5<-fit_UMW(Useful$EMBORCAÇÃO)
AjSeMw6<-fit_UMW(Useful$FURNAS)
AjSeMw7<-fit_UMW(Useful$ITUMBIARA)
AjSeMw8<-fit_UMW(Useful$NOVA.PONTE)
AjSeMw9<-fit_UMW(Useful$SÃO.SIMÃO)
AjSeMw10<-fit_UMW(Useful$S.DO.FACÃO)

### Northeast/North ----

AjNe1<-fit_UMW(Useful$IRAPE)
AjNe2<-fit_UMW(Useful$ITAPARICA)
AjNe3<-fit_UMW(Useful$SOBRADINHO)
AjNe4<-fit_UMW(Useful$TRÊS.MARIAS)
AjN5<-fit_UMW(Useful$SERRA.DA.MESA)
AjN6<-fit_UMW(Useful$TUCURUI)

## Reading Skills --------------------------------------------------------------

data3<-read.table("data_set/reading_skills.txt", header = TRUE)
data3$intera2<-data3$dislex*(data3$QI2)

y<-data3$resp
X<-data3[,-1]
X1<-X[,c(3,5)]

fit<-fit_RQUMW(y = y,X = X1,tau = 0.5,link = "logit",intercept = T,printmodel = T)
# summary_rqumw(fit)

### Test cov -------------------------------------------------------------------

cores <- parallel::detectCores()
cl <- parallel::makeCluster(cores-1) # not to overload your computer
doSNOW::registerDoSNOW(cl)
# Secelct ---
library(foreach)
models<-suppressWarnings(select_cov(writing = "dislex",ncovn=c(1,2),tau = 0.5,
                                    X=X,y=y,criterion = AIC,link = "logit",
                                    intercept = T,save = F))
models
# stop parallel ---
foreach::registerDoSEQ()
parallel::stopCluster(cl)


