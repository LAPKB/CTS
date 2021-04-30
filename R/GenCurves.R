# Script to generate curves for clinical trial simulator examples.
# MN and JB. Modified 4/20/21.



#' @export
gen_Sample <- function(){
  require(tidyverse)
  require(Pmetrics)
  #simulate N time concentration curves

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
  getAUC <- function(x){
    auc <- makeAUC(x,conc~time)
    return(auc)
  }
  
  #get Cmax
  getCmax <- function(x){
    cmax <- tapply(x$conc,x$id,max)
    return(cmax)
  }
  
  #add IOV
  addIOV <- function(x,iov){
    x$conc <- unlist(tapply(x$conc,x$id,
                            function(h) h * exp(rnorm(1,mean=0,sd=iov))))
    return(x)
  }

  testBE <- function(ref,comp,paired){
    t_res <- t.test(log10(ref),log10(comp),paired=paired,conf.level=0.9)
    return(10**t_res$conf.int)
    
  }
  
  
  # MAKE PROFILES -----------------------------------------------------------
  set.seed(1)
  
  #reference profiles
  ref <- simCurves(1000,meanPar = c(log(.2),log(5)), 
                   sdPar = c(log(2),log(2)),dose=100)
  
  # #generic 1 profiles (truly BE)
  # gen1 <- simCurves(1000,meanPar = c(log(.22),log(5.1)), 
  #                   sdPar = c(log(2.2),log(2.2)),dose=100)
  
  #generic 2 profiles (truly NOT BE)
  gen1 <- simCurves(1000,meanPar = c(log(.4),log(5.2)),
                    sdPar = c(log(1.2),log(2.2)),dose=100)
  
  iov.sd=.3 
  ref_iov <- addIOV(ref,iov.sd)
  gen1_iov <- addIOV(gen1,iov.sd)
  gen2_iov <- addIOV(gen2,iov.sd)
  
  
  
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
  
  
  #make sample data
  
  #Ref is drug A
  # auc_ref$seq[1:500] <- 1 # ABAB
  # auc_ref$seq[501:1000] <- 2 #BABA
  # auc_ref$per[1:500] <- 1
  # auc_ref$per[501:1000] <- 2 
  auc_ref$per <- 1
  auc_ref$drug <- "A"
  
  #Gen1 is drug B
  # auc_gen1$seq[1:500] <- 1 # ABAB
  # auc_gen1$seq[501:1000] <- 2 #BABA
  # auc_gen1$per[1:500] <- 2
  # auc_gen1$per[501:1000] <- 1 
  auc_gen1$per <- 1
  auc_gen1$drug <- "B"
  
  #Ref is drug A
  # auc_ref_iov$seq[1:500] <- 1 # ABAB
  # auc_ref_iov$seq[501:1000] <- 2 #BABA
  # auc_ref_iov$per[1:500] <- 3
  # auc_ref_iov$per[501:1000] <- 4 
  auc_ref_iov$per <- 2
  auc_ref_iov$drug <- "A"
  
  #Gen1 is drug B
  # auc_gen1_iov$seq[1:500] <- 1 # ABAB
  # auc_gen1_iov$seq[501:1000] <- 2 #BABA
  # auc_gen1_iov$per[1:500] <- 4
  # auc_gen1_iov$per[501:1000] <- 3 
  auc_gen1_iov$per <- 2
  auc_gen1_iov$drug <- "B"
  
  sampData <- as_tibble(rbind(auc_ref,auc_ref_iov,auc_gen1,auc_gen1_iov)) %>%
    mutate(cmax = c(cmax_ref,cmax_ref_iov,cmax_gen1,cmax_gen1_iov)) %>%
    mutate(across(.cols=c("id","per","drug"),factor)) %>%
    select(id, auc=tau, cmax, per, drug)
  
  return(sampData)
  
}



s1 <- gen_Sample()
