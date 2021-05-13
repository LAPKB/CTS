#' Generate simulated sample data
#' 
#' Creates concentration time curves for clinical trial simulator examples from a 
#' one-compartment model with instantaneous bolus input.
#' 
#' @param N number of profiles per condition
#' @param meanPar List of length two. Each item is a vector with the log10(mean) 
#' for Ke and V.  The first item is the reference compound; the second item is the
#' generic.
#' @param sdPar List of length two. Each item is a vector with the log10(sd) 
#' for Ke and V.  The first item is the reference compound; the second item is the
#' generic.
#' @param dose Simulated dose for both reference and generic compounds
#' @param iov.sd SD to apply for random intra-occasion variability in parameter values.
#' @param seed Seed for random number generator.
#' @return A tibble with columns ID, time, conc, per, drug
#' @author Michael Neely
#' @export

gen_Sample <- function(N=1000,meanPar = list(c(log(.2),log(5)), c(log(.22),log(4.5))),
                       sdPar = list(c(log(.5),log(.5)), c(log(3),log(3))),
                       dose = 100, iov.sd = 0.5, seed = -17){

  #simulate N time concentration curves
  
  simCurves <- function(N, thisMean, thisSD, dose){
    #parameters should be in log scale
    KE_mu<- thisMean[1]
    V_mu <- thisMean[2]
    KE_sd <- abs(thisSD[1])
    V_sd <- abs(thisSD[2])
    times <- seq(0,24,1)
    set.seed(seed)
    KE <- rlnorm(N,KE_mu,KE_sd)
    V <- rlnorm(N,V_mu,V_sd)
    conc <- data.frame(sapply(times, function(t) dose/V * exp(-KE*t)))
    conc2 <- pivot_longer(conc,cols=everything())
    names(conc2) <- c("id","conc")
    conc2$id <- rep(1:N,each=length(times))
    conc2$time <- rep(times,N)
    conc2 <- conc2[,c("id","time","conc")]
    return(data.frame(conc2))
  }
  
  #add IOV
  addIOV <- function(x,iov){
    x$conc <- unlist(tapply(x$conc,x$id,
                            function(h) h * exp(rnorm(1,mean=0,sd=iov))))
    return(x)
  }
  
  
  
  # MAKE PROFILES -----------------------------------------------------------
 
  #reference profiles
  ref <- simCurves(N,meanPar[[1]],sdPar[[1]],dose)
  ref$per <- 1
  ref$drug <- "R"
  
  #generic profiles 
  gen1 <- simCurves(N,meanPar[[2]],sdPar[[2]],dose)
  gen1$per <- 1
  gen1$drug <- "T1"
  
  #iov
  ref_iov <- addIOV(ref,iov.sd)
  ref_iov$per <- 2
  ref_iov$drug <- "R"
  gen1_iov <- addIOV(gen1,iov.sd)
  gen1$per <- 2
  gen1$drug <- "T1"
  
  sampData <- as_tibble(rbind(ref,ref_iov,gen1,gen1_iov)) %>%
    mutate(across(.cols=c("id","per","drug"),factor)) 
  
  return(sampData)

}
  
 
 


