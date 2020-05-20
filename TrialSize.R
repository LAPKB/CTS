



ratiomean <- vector("numeric")
lower <- vector("numeric")
upper <- vector("numeric")
ssize <- 5:40
k <- 0
for(i in ssize){
  ratios <- list()
  k <- k+1
  for(j in 1:i){
    ratios[[j]] <- rnorm(1000,mean=runif(1,-.1,.25),sd=0.1)
  }
  temp <- sapply(ratios,mean)
  lower[k] <- 10**t.test(temp,mu=0,conf.level=0.9)$conf.int[1]
  upper[k] <- 10**t.test(temp,mu=0,conf.level=0.9)$conf.int[2]
  
  ratiomean[k] <- mean(temp)
  
}

plot(x=ssize,y=upper,type="o",ylim=c(min(c(0.8,lower)),max(c(1.25,upper))),xlab="Sample Size per Group",ylab="GMR 90% CI",cex.lab=1.2)
lines(x=ssize,y=lower,type="o",lty=2)
abline(h=c(0.8,1.25),col=c("red","blue"))




#simulate N time concentration curves
library(tidyr)
library(ggplot2)
library(Pmetrics)

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

#reference profiles
ref <- simCurves(1000,meanPar = c(log(.2),log(5)), sdPar = c(log(2),log(2)),dose=100)
#generic 1 profiles (truly BE)
gen1 <- simCurves(1000,meanPar = c(log(.22),log(5.5)), sdPar = c(log(2),log(2)),dose=100)
#generic 2 profiles (truly NOT BE)
gen2 <- simCurves(1000,meanPar = c(log(.1),log(4)), sdPar = c(log(2),log(2)),dose=100)

#plot curves
p <- ggplot(ref,aes(x=time,y=conc,group=id)) + geom_line()
p
p2 <- ggplot(gen1,aes(x=time,y=conc,group=id)) + geom_line()
p2
p3 <- ggplot(gen2,aes(x=time,y=conc,group=id)) + geom_line()
p3

#get AUC
getAUC <- function(x){
  auc <- makeAUC(x,conc~time)
  return(auc)
}

auc_ref <- getAUC(ref)
auc_gen1 <- getAUC(gen1)
auc_gen2 <- getAUC(gen2)

#overall check
#CI should be [-0.2, 0.2]
t.test(log10(auc_ref$tau),log10(auc_gen1$tau),paired=T,conf.level=0.9) #yes, BE
t.test(log10(auc_ref$tau),log10(auc_gen2$tau),paired=T,conf.level=0.9) #no, not BE


