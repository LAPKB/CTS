# Script for simulation of parallel design using NPOD-generated curves.
# Jay Bartroff, last modified 11.Mar.20.

library(tidyr)
library(ggplot2)
library(Pmetrics)


# FUNCTIONS ---------------------------------------------------------------

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

#'@export
testBE <- function(ref,comp,paired){
  t_res <- t.test(log10(ref),log10(comp),paired=paired,conf.level=0.9)
  return(10**t_res$conf.int)
  
}

#TODO: we need to get this into functions

# Read in NPOD data files
ref=read.csv(file="data/ref_data.csv", sep=",")
gen1=read.csv(file="data/test_data.csv", sep=",")

# strip occasion, adjust column names
ref=ref[,-4]; names(ref) = c("id", "time", "conc") 
gen1=gen1[,-4]; names(gen1) = c("id", "time", "conc")
  

#plot curves as sanity check
ggplot(ref,aes(x=time,y=conc,group=id)) + geom_line()
ggplot(gen1,aes(x=time,y=conc,group=id)) + geom_line()

# EXTRACT AUC AND CMAX ----------------------------------------------------

auc_ref <- getAUC(ref)
auc_gen1 <- getAUC(gen1)

cmax_ref <- getCmax(ref)
cmax_gen1 <- getCmax(gen1)


#overall checks
testBE(auc_ref$tau,auc_gen1$tau,paired=TRUE) 
testBE(cmax_ref,cmax_gen1,paired=T) 

#####################
#  Parallel design.  

set.seed(110321)
n.curves=1000
n.trials=100  # per sample size
#n.subj=seq(from=5,to=40, by=5) # no. subjects per sequence
n.subj=seq(from=50,to=250, by=10)
no.subj=length(n.subj)
pBE1=rep(0,times=no.subj)
up.s1=rep(0,times=no.subj); up.ss1=rep(0,times=no.subj)
dn.s1=rep(0,times=no.subj); dn.ss1=rep(0,times=no.subj)
pBE2=rep(0,times=no.subj)
up.s2=rep(0,times=no.subj); up.ss2=rep(0,times=no.subj)
dn.s2=rep(0,times=no.subj); dn.ss2=rep(0,times=no.subj)
for (i in 1:no.subj) {
  for (b in 1:n.trials) {
    subjs=sample.int(n=n.curves, size=2*n.subj[i])
    seq1=head(subjs,n.subj[i]); seq2=tail(subjs,n.subj[i])
    ci=testBE(auc_ref$tau[seq1],auc_gen1$tau[seq2],paired=F)
    up.s1[i]=up.s1[i]+ci[2]; up.ss1[i]=up.ss1[i]+(ci[2])^2
    dn.s1[i]=dn.s1[i]+ci[1]; dn.ss1[i]=dn.ss1[i]+(ci[1])^2
    if ((ci[1]>=.8) & (ci[2]<=1.25)) {
      pBE1[i]=pBE1[i]+1
    }
  }
}

pBE1=pBE1/n.trials
up.s1=up.s1/n.trials; up.ss1=up.ss1/n.trials
dn.s1=dn.s1/n.trials; dn.ss1=dn.ss1/n.trials



# Plots
#par(mfrow=c(2,2))
plot(x=n.subj, y=pBE1, main="Prob of BE: Ref/Test BE", ylim=c(0,1), type="l", ylab="Probability", xlab="Sample size per sequence")
p1.sd=sqrt(pBE1*(1-pBE1)/n.trials)
arrows(x0=n.subj, y0=pBE1-p1.sd, x1=n.subj, y1=pBE1+p1.sd, code=3, angle=90, length=0.02, col="blue")

# plot(x=n.subj, y=up.s1, main="Ref/Test BE", ylim=c(.5,2),xlab="Sample size per sequence")
# points(x=n.subj, y=dn.s1)
plot(1, type="n",main="Avg Effect Ratio: Ref/Test BE", xlab="Sample size per sequence", ylab="", xlim=c(min(n.subj), max(n.subj)), ylim=c(.6,1.75))
abline(h=.8, col="red", lty="dashed"); abline(h=1.25, col="red", lty="dashed")
arrows(x0=n.subj, y0=dn.s1, x1=n.subj, y1=up.s1, code=3, angle=90, length=0.02)
