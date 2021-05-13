

library(dplyr)
library(tidyr)
library(ggplot2)
library(Pmetrics)
library(openxlsx)
library(readr)

# FUNCTIONS ---------------------------------------------------------------

#simulate N time concentration curves
#'@export
simCurves <- function(N,meanPar,sdPar,dose){
  #parameters should be in log scale
  KE_mu<- meanPar[1]
  V_mu <- meanPar[2]
  KE_sd <- sdPar[1]
  V_sd <- sdPar[2]
  times <- seq(0,24,1)
  set.seed(-1)
  KE <- rlnorm(N,KE_mu,KE_sd)
  V <- rlnorm(N,V_mu,V_sd)
  conc <- data.frame(sapply(times, function(t) dose/V * exp(-KE*t)))
  conc2 <- pivot_longer(conc,cols=everything())
  names(conc2) <- c("id","conc")
  conc2$id <- rep(1:1000,each=length(times))
  conc2$time <- rep(times,N)
  conc2 <- conc2[,c("id","time","conc")]
  return(data.frame(conc2))
}

#get AUC
#'@export
getAUC <- function(x){
  auc <- makeAUC(x,conc~time)
  return(auc)
}

#get Cmax
#'@export
getCmax <- function(x){
  cmax <- tapply(x$conc,x$id,max)
  return(cmax)
}

#add IOV
#'@export
addIOV <- function(x,iov){
  x$conc <- unlist(tapply(x$conc,x$id,function(h) h * exp(rnorm(1,mean=0,sd=iov))))
  return(x)
}

#'@export
testBE <- function(ref,comp,paired){
  t_res <- t.test(log10(ref),log10(comp),paired=paired,conf.level=0.9)
  return(10**t_res$conf.int)
  
}

#TODO: get this into functions

# MAKE PROFILES -----------------------------------------------------------

#reference profiles
ref <- simCurves(1000,meanPar = c(log(.2),log(5)), sdPar = c(log(2),log(2)),dose=100)
#generic 1 profiles (truly BE)
gen1 <- simCurves(1000,meanPar = c(log(.22),log(5.1)), sdPar = c(log(2.2),log(2.2)),dose=100)
#generic 2 profiles (truly NOT BE)
gen2 <- simCurves(1000,meanPar = c(log(.22),log(6)), sdPar = c(log(7),log(7)),dose=100)

ref_iov <- addIOV(ref,0.3)
gen1_iov <- addIOV(gen1,0.3)
gen2_iov <- addIOV(gen2,0.3)


#plot curves as sanity check
ggplot(ref,aes(x=time,y=conc,group=id)) + geom_line()
ggplot(gen1,aes(x=time,y=conc,group=id)) + geom_line()
ggplot(gen2,aes(x=time,y=conc,group=id)) + geom_line()


ggplot(ref_iov,aes(x=time,y=conc,group=id)) + geom_line()
ggplot(gen1_iov,aes(x=time,y=conc,group=id)) + geom_line()
ggplot(gen2_iov,aes(x=time,y=conc,group=id)) + geom_line()



# EXTRACT AUC AND CMAX ----------------------------------------------------

auc_ref <- getAUC(ref)
auc_gen1 <- getAUC(gen1)
auc_gen2 <- getAUC(gen2)
auc_ref_iov <- getAUC(ref_iov)
auc_gen1_iov <- getAUC(gen1_iov)
auc_gen2_iov <- getAUC(gen2_iov)

cmax_ref <- getCmax(ref)
cmax_gen1 <- getCmax(gen1)
cmax_gen2 <- getCmax(gen2)
cmax_ref_iov <- getCmax(ref_iov)
cmax_gen1_iov <- getCmax(gen1_iov)
cmax_gen2_iov <- getCmax(gen2_iov)

#overall checks
#CI should be [0.8, 1.25]

testBE(auc_ref$tau,auc_gen1$tau,paired=T) #yes, BE
testBE(auc_ref$tau,auc_gen2$tau,paired=T) #no, not BE

testBE(cmax_ref,cmax_gen1,paired=T) #yes, BE
testBE(cmax_ref,cmax_gen2,paired=T) #no, not BE


#try with buproprion data
sr <- read.xlsx("data/Bup150SR.xlsx",na.strings=".") %>% 
  pivot_longer(cols=-time,names_to = "id", values_to = "conc") %>%
  arrange(id,time) %>% 
  select(id,time,conc) %>% 
  filter(!is.na(conc))

xl<- read.xlsx("data/Bup150XL.xlsx",na.strings=".") %>% 
  pivot_longer(cols=-time,names_to = "id", values_to = "conc") %>%
  arrange(id,time) %>% 
  select(id,time,conc) %>% 
  filter(!is.na(conc))




srAUC <- getAUC(data.frame(sr))
xlAUC <- getAUC(data.frame(xl))

srAUC <- srAUC[srAUC$id %in% xlAUC$id,]
xlAUC <- xlAUC[xlAUC$id %in% srAUC$id,]


testBE(xlAUC$tau,srAUC$tau,paired=T) #yes, BE


#2x2 parallel
#randomly select N=5 profiles from each group (ref vs gen1 OR ref vs gen2)
#this simulates different individuals
#testBE for AUC and Cmax
#repeat 100 times (100 clinical trials)
#store 1) percent of 100 trials with 90%CI GMR [0.8, 1.25] and 2) mean upper and lower 90% CI
#repeat for N = 10, 15, 20, 25, 30, 35, 40

#2x2 crossover
#as for parallel but match selection in each group, i.e. select 5 random curves from ref
#and same curves from comparator
#this simulates same individual getting each drug

#2x2 crossover with replication (partial or full)
#as for crossover but use IOV curves for second dosing event in subjects
#appropriate statistics to account for replication


