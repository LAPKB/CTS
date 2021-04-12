# Script for simulation of 2x2x4 crossover design.  Assumes existence in the workspace of auc_ref, auc_gen1, auc_gen2, auc_ref_iov, auc_gen1_iov, auc_gen2_iov as outputs of GenCurves.R, or some other source.
# Jay Bartroff, last modified 7.Apr.20.

library(nlme)

ctrl=lmeControl(opt='optim', msMaxIter=1000)
ci.lcut=.8; ci.ucut=1.25 # lower and upper cutoffs for CI

set.seed(3)
n.curves=1000
n.trials=50  # per sample size
n.subj=seq(from=5,to=25, by=5) # no. subjects per sequence, the max needs to be <= n.curves/2
#n.subj=50 # a single SS for debugging
n.SS=length(n.subj)
p.dropout=.3 # currently, this fraction of lines of the data frame are randomly dropped

pBE1=rep(0,times=n.SS)
rat.s1=rep(0,times=n.SS) # estimate of AUC ratio; so SoS because recording CI too
up.s1=rep(0,times=n.SS)#; up.ss1=rep(0,times=n.SS)
dn.s1=rep(0,times=n.SS)#; dn.ss1=rep(0,times=n.SS)
dd.eff.s1=rep(0,times=n.SS); dd.eff.ss1=rep(0,times=n.SS) #direct drug effect "test-ref" (actually diff of logs)
intra.CV.s1=rep(0,times=n.SS); intra.CV.ss1=rep(0,times=n.SS) # Intra-subject CV

pBE2=rep(0,times=n.SS)
rat.s2=rep(0,times=n.SS) # estimate of AUC ratio; so SoS because recording CI too
up.s2=rep(0,times=n.SS)#; up.ss2=rep(0,times=n.SS)
dn.s2=rep(0,times=n.SS)#; dn.ss2=rep(0,times=n.SS)
dd.eff.s2=rep(0,times=n.SS); dd.eff.ss2=rep(0,times=n.SS)
intra.CV.s2=rep(0,times=n.SS); intra.CV.ss2=rep(0,times=n.SS)

# For building column vectors to make up the data frame for each simulated trial, my convention is to move down rows in the design layout (e.g., Chow & Liu, p. 40), then over columns.  These columns will be seq, prd (period), drug, and AUC (or other observation, e.g., Cmax). I use the BEAR convention (see http://pkpd.kmu.edu.tw/bear/) of seq=1 is where the ref is taken first and  drug=1 means reference (drug=2 means test, e.g., generic).

for (i in 1:n.SS) {
  cat('SS ',i,' of ',n.SS,sep='','\n')
  n.rows.df=8*n.subj[i] # no. rows in d.f. (or vectors if not using d.f.)
  
  # make d.f. columns that are static over the simulated trials for this SS
  seq=as.factor(rep.int(c(rep.int(1,n.subj[i]),rep.int(2,n.subj[i])),4))
  prd=as.factor(c(rep.int(1,2*n.subj[i]), rep.int(2,2*n.subj[i]), rep.int(3,2*n.subj[i]), rep.int(4,2*n.subj[i])))
  drug=as.factor(c(rep.int(1,n.subj[i]), rep.int(2,2*n.subj[i]), rep.int(1,2*n.subj[i]), rep.int(2,2*n.subj[i]), rep.int(1,n.subj[i])))
  
  for (b in 1:n.trials) {
    cat('   Trial ',b,' of ',n.trials,sep='','\n')
    
    #*************
    # Comparison 1
    subjs=sample.int(n=n.curves, size=2*n.subj[i])
    seq1=head(subjs,n.subj[i]); seq2=tail(subjs,n.subj[i])
  
    # make subject and observation columns for d.f., which change each trial
    subj=as.factor(rep.int(c(seq1,seq2),4))
    AUC.log1=log(c(auc_ref$tau[seq1], auc_gen1$tau[seq2], auc_gen1$tau[seq1], auc_ref$tau[seq2], auc_ref_iov$tau[seq1], auc_gen1_iov$tau[seq2], auc_gen1_iov$tau[seq1], auc_ref_iov$tau[seq2]))
    #this.df=data.frame(subj,seq,prd,drug,AUC.log)
    
    # mixed model
    rows2incl=1:n.rows.df
    if (p.dropout>0) {
      rows2incl=rows2incl[-sample.int(n=n.rows.df, size=round(p.dropout*n.rows.df))]
    }
    #print(rows2incl)
      
    target.lme1=lme(AUC.log1 ~ seq +  prd + drug, random=~drug - 1|subj, control=ctrl, weights=varIdent(form = ~ 1 | drug), method="REML", na.action=na.exclude, subset=rows2incl) 
    
    # record stuff
    dd.eff1=summary(target.lme1)[20][[1]][6,1]  # direct drug effect.  BTW, index 20 of summary() is the coefficient t-table, and the 6th row is for the drug2 effect: Value, Std.Error, DF, t-value, p-value
    df1=summary(target.lme1)[20][[1]][6,3]
    se1=summary(target.lme1)[20][[1]][6,2]
    rat1=exp(dd.eff1)
    deltaCI1=qt(.05,df1,lower.tail=FALSE)*se1
    ci.l1=rat1*exp(-deltaCI1)
    ci.u1=rat1*exp(deltaCI1)
    if ((ci.l1>=ci.lcut) & (ci.u1<=ci.ucut)) pBE1[i]=pBE1[i]+1
    rat.s1[i]=rat.s1[i]+rat1
    dn.s1[i]=dn.s1[i]+ci.l1
    up.s1[i]=up.s1[i]+ci.u1
    dd.eff.s1[i]=dd.eff.s1[i]+dd.eff1
    dd.eff.ss1[i]=dd.eff.ss1[i]+dd.eff1^2
    mse1=2*(deltaCI1/((sqrt(2/n.subj[i])*qt(.05,2*n.subj[i]-2,lower.tail=FALSE))))^2
    cv.intra1=sqrt(exp(mse1)-1)
    # Bear code for this is:
    # MSE <- 2*(delta_CI/((sqrt(1/L1+1/L2)*qt(0.95,L1+L2-2))))^2
    # CVintra <- 100*sqrt(exp(MSE)-1)
    intra.CV.s1[i]=intra.CV.s1[i]+cv.intra1
    intra.CV.ss1[i]=intra.CV.ss1[i]+cv.intra1^2
    
    #*************
    # Comparison 2
    subjs=sample.int(n=n.curves, size=2*n.subj[i])
    seq1=head(subjs,n.subj[i]); seq2=tail(subjs,n.subj[i])
    
    # make subject and observation columns for d.f., which change each trial
    subj=as.factor(rep.int(c(seq1,seq2),4))
    AUC.log2=log(c(auc_ref$tau[seq1], auc_gen2$tau[seq2], auc_gen2$tau[seq1], auc_ref$tau[seq2], auc_ref_iov$tau[seq1], auc_gen2_iov$tau[seq2], auc_gen2_iov$tau[seq1], auc_ref_iov$tau[seq2]))
    
    # mixed model
    rows2incl=1:n.rows.df
    if (p.dropout>0) {
      rows2incl=rows2incl[-sample.int(n=n.rows.df, size=round(p.dropout*n.rows.df))]
    }
    #print(rows2incl)
    
    target.lme2=lme(AUC.log2 ~ seq +  prd + drug, random=~drug - 1|subj, control=ctrl, weights=varIdent(form = ~ 1 | drug), method="REML", na.action=na.exclude, subset=rows2incl) 
    
    # record stuff
    dd.eff2=summary(target.lme2)[20][[1]][6,1]  # direct drug effect.  BTW, index 20 of summary() is the coefficient t-table, and the 6th row is for the drug2 effect: Value, Std.Error, DF, t-value, p-value
    df2=summary(target.lme2)[20][[1]][6,3]
    se2=summary(target.lme2)[20][[1]][6,2]
    rat2=exp(dd.eff2)
    deltaCI2=qt(.05,df2,lower.tail=FALSE)*se2
    ci.l2=rat2*exp(-deltaCI2)
    ci.u2=rat2*exp(deltaCI2)
    if ((ci.l2>=ci.lcut) & (ci.u2<=ci.ucut)) pBE2[i]=pBE2[i]+1
    rat.s2[i]=rat.s2[i]+rat2
    dn.s2[i]=dn.s2[i]+ci.l2
    up.s2[i]=up.s2[i]+ci.u2
    dd.eff.s2[i]=dd.eff.s2[i]+dd.eff2
    dd.eff.ss2[i]=dd.eff.ss2[i]+dd.eff2^2
    mse2=2*(deltaCI2/((sqrt(2/n.subj[i])*qt(.05,2*n.subj[i]-2,lower.tail=FALSE))))^2
    cv.intra2=sqrt(exp(mse2)-1)
    # Bear code for this is:
    # MSE <- 2*(delta_CI/((sqrt(1/L1+1/L2)*qt(0.95,L1+L2-2))))^2
    # CVintra <- 100*sqrt(exp(MSE)-1)
    intra.CV.s2[i]=intra.CV.s2[i]+cv.intra2
    intra.CV.ss2[i]=intra.CV.ss2[i]+cv.intra2^2
    
  }
}

pBE1=pBE1/n.trials
rat.s1=rat.s1/n.trials
up.s1=up.s1/n.trials
dn.s1=dn.s1/n.trials
dd.eff.s1=dd.eff.s1/n.trials; dd.eff.ss1=dd.eff.ss1/n.trials
intra.CV.s1=intra.CV.s1/n.trials; intra.CV.ss1=intra.CV.ss1/n.trials

pBE2=pBE2/n.trials
rat.s2=rat.s2/n.trials
up.s2=up.s2/n.trials
dn.s2=dn.s2/n.trials
dd.eff.s2=dd.eff.s2/n.trials; dd.eff.ss2=dd.eff.ss2/n.trials
intra.CV.s2=intra.CV.s2/n.trials; intra.CV.ss2=intra.CV.ss2/n.trials

## Make lots of various plots

# Prob(BE)
par(mfrow=c(2,1))

plot(x=n.subj, y=pBE1, main="Prob of BE: Truly BE", ylim=c(0,1), type="l", ylab="Probability", xlab="Sample size per sequence")
p1.sd=sqrt(pBE1*(1-pBE1)/n.trials)
arrows(x0=n.subj, y0=pBE1-p1.sd, x1=n.subj, y1=pBE1+p1.sd, code=3, angle=90, length=0.02, col="blue")

plot(x=n.subj, y=pBE2, main="Prob of BE: Not BE", ylim=c(0,1), type="l", ylab="Probability", xlab="Sample size per sequence")
p2.sd=sqrt(pBE2*(1-pBE2)/n.trials)
arrows(x0=n.subj, y0=pBE2-p2.sd, x1=n.subj, y1=pBE2+p2.sd, code=3, angle=90, length=0.02, col="blue")

# Ratio and CI
par(mfrow=c(2,1))

plot(1, type="n",main="Avg Effect Ratio: Truly BE", xlab="Sample size per sequence", ylab="", xlim=c(min(n.subj), max(n.subj)), ylim=c(.5,1.5))
abline(h=.8, col="red", lty="dashed"); abline(h=1.25, col="red", lty="dashed")
arrows(x0=n.subj, y0=dn.s1, x1=n.subj, y1=up.s1, code=3, angle=90, length=0.02, col="blue")
points(x=n.subj, y=rat.s1, pch=4)

plot(1, type="n",main="Avg Effect Ratio: Not BE", xlab="Sample size per sequence", ylab="", xlim=c(min(n.subj), max(n.subj)), ylim=c(.1,2))
abline(h=.8, col="red", lty="dashed"); abline(h=1.25, col="red", lty="dashed")
arrows(x0=n.subj, y0=dn.s2, x1=n.subj, y1=up.s2, code=3, angle=90, length=0.02, col="blue")
points(x=n.subj, y=rat.s2, pch=4)

# Direct drug effect
# both on same
plot(1, type="n",main="Direct Drug Effect", xlab="Sample size per sequence", ylab="log(test/ref)", xlim=c(min(n.subj), max(n.subj)), ylim=c(-.25,-.08))#ylim=c(-.13,-.08))
arrows(x0=n.subj, x1=n.subj, y0=dd.eff.s1-sqrt((dd.eff.ss1-dd.eff.s1^2)/n.trials), y1=dd.eff.s1+sqrt((dd.eff.ss1-dd.eff.s1^2)/n.trials), code=3, angle=90, length=0.02, col="blue")
points(x=n.subj, y=dd.eff.s1, pch=4)
arrows(x0=n.subj, x1=n.subj, y0=dd.eff.s2-sqrt((dd.eff.ss2-dd.eff.s2^2)/n.trials), y1=dd.eff.s2+sqrt((dd.eff.ss2-dd.eff.s2^2)/n.trials), code=3, angle=90, length=0.02, col="red")
points(x=n.subj, y=dd.eff.s2, pch=4)
legend(x=20,y=-.15, c("Truly BE","Not BE"), cex=1,col=c("blue","red"), lty=1)
#####

# separate
plot(1, type="n",main="Direct Drug Effect: Truly BE", xlab="Sample size per sequence", ylab="log(test/ref)", xlim=c(min(n.subj), max(n.subj)), ylim=c(-.37,-.08))#ylim=c(-.13,-.08))
arrows(x0=n.subj, x1=n.subj, y0=dd.eff.s1-sqrt((dd.eff.ss1-dd.eff.s1^2)/n.trials), y1=dd.eff.s1+sqrt((dd.eff.ss1-dd.eff.s1^2)/n.trials), code=3, angle=90, length=0.02, col="blue")
points(x=n.subj, y=dd.eff.s1, pch=4)

plot(1, type="n",main="Direct Drug Effect: Not BE", xlab="Sample size per sequence", ylab="log(test/ref)", xlim=c(min(n.subj), max(n.subj)), ylim=c(-.37,-.08)) #ylim=c(-.37,-.2))
arrows(x0=n.subj, x1=n.subj, y0=dd.eff.s2-sqrt((dd.eff.ss2-dd.eff.s2^2)/n.trials), y1=dd.eff.s2+sqrt((dd.eff.ss2-dd.eff.s2^2)/n.trials), code=3, angle=90, length=0.02, col="blue")
points(x=n.subj, y=dd.eff.s2, pch=4)

# Intra subject CV
# both on same
plot(1, type="n",main="Intra Subject CV for log(test/ref)", xlab="Sample size per sequence", ylab="", xlim=c(min(n.subj), max(n.subj)), ylim=c(.15,.185))
arrows(x0=n.subj, x1=n.subj, y0=intra.CV.s1-sqrt((intra.CV.ss1-intra.CV.s1^2)/n.trials), y1=intra.CV.s1+sqrt((intra.CV.ss1-intra.CV.s1^2)/n.trials), code=3, angle=90, length=0.02, col="blue")
points(x=n.subj, y=intra.CV.s1, pch=4)
arrows(x0=n.subj, x1=n.subj, y0=intra.CV.s2-sqrt((intra.CV.ss2-intra.CV.s2^2)/n.trials), y1=intra.CV.s2+sqrt((intra.CV.ss2-intra.CV.s2^2)/n.trials), code=3, angle=90, length=0.02, col="red")
points(x=n.subj, y=intra.CV.s2, pch=4)
legend(x=20,y=.16, c("Truly BE","Not BE"), cex=1,col=c("blue","red"), lty=1)
#####


# Function for 5%, 50%, and 95% time-concentration curves, computed directly from the database.
TC.percs = function(tc.table, percs4TC = c(.05, .5, .95), ylim.plot=NULL, main.plot="Percentiles of Time-Concentration Curves", cols.plot=rainbow(length(percs4TC)), leg.loc="topright", lwd.plot=4, alpha.f.plot=.075) {
  # tc.table = some data frame like ref, etc., above
  # percs4TC = which percentiles to keep track of
  #  J.B. 18.Nov.20; last modified same
  
  n.percsTC = length(percs4TC)
  time.pts = unique(tc.table$time[tc.table$id == 1])
  n.time.pts=length(time.pts)
  percs.dat.TC = matrix(nrow = n.time.pts, ncol = n.percsTC)  # percs.dat.TC[i,j]=j-th percentile we're keeping track of, at i-th time point
  for (i in 1:n.time.pts) {
    percs.dat.TC[i,]=quantile(x=tc.table$conc[tc.table$time==tc.table$time[i]], probs=percs4TC, names=FALSE)
  }
  
  # plot
  if (is.null(ylim.plot)) ylim.plot=c(min(percs.dat.TC), max(percs.dat.TC))
  plot(x=time.pts, y=percs.dat.TC[,1], ylim=ylim.plot, xlab="Time", ylab = "Concentration", main=main.plot, col=cols.plot[1], lty=1, type="l", lwd=lwd.plot)
  for (j in 2:n.percsTC) {
    lines(x=time.pts, y=percs.dat.TC[,j], col=cols.plot[j], lwd=lwd.plot)
  }
  legend(x=leg.loc, legend=paste(100*percs4TC,"%", sep=''), fill=cols.plot)
  
  # plot individual curves underneath
  indivs=unique(tc.table$id)
  for (indiv in indivs) {
    rws=(tc.table$id==indiv)
    lines(x=tc.table$time[rws], y=tc.table$conc[rws], col = adjustcolor('black', alpha.f=alpha.f.plot))
  }
  
  return(percs.dat.TC)
}

TC.percs(tc.table=ref, main.plot="Reference (Database)", lwd.plot=5, cols.plot=c("red", "black", "blue"), ylim.plot=c(0,80))
TC.percs(tc.table=gen1, main.plot="Generic 1 - Truly BE (Database)", lwd.plot=5, cols.plot=c("red", "black", "blue"), ylim.plot=c(0,80))
TC.percs(tc.table=gen2, main.plot="Generic 2 - Not BE (Database)", lwd.plot=5, cols.plot=c("red", "black", "blue"), ylim.plot=c(0,80))

