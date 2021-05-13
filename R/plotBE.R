#' Plot BE results
#' 
#' Generates two types of plots: ratio and probability of BE for given sample size.
#' Ratio plots are the geometric mean by sample size of either AUC or Cmax and 
#' confidence interval, as calculated by \code{\link{extractBE}}. 
#' Probability plots are the probability of demonstrating BE by sample size.
#' 
#' @param x Output of \code{\link{extractBE}} to plot.
#' @param ci_lo Lower confidence interval reference line for geometric mean AUC or Cmax to be
#' considered BE, default 0.8 per FDA guidance.
#' @param ci_up Upper confidence interval reference line for geometric mean AUC or Cmax to be
#' considered BE, default 1.25 per FDA guidance.
#' @param ylab Label for y-axis to override defaults.
#' @return A plot.
#' @author Michael Neely
#' @export

plotBE <- function(x,type, ci_lo = 0.8, ci_up = 1.25, ylab){
  if(all(is.na(x))) {
    cat("Replicated design not simulated.\n")
    return(invisible())
  }
  if(type == "ratio"){
    if(missing(ylab)) ylab <- "Ratio"
    p <- ggplot(x,aes(x=n_subj,y=de_ratio)) + geom_line(lwd=1.2) + 
      geom_errorbar(aes(ymin=ci_lo, ymax=ci_up, width=0.3), col = "dodgerblue", alpha = 0.5) + #ylim(c(0.5,1.3)) +
      geom_hline(yintercept = c(ci_lo, ci_up), lty = 2, col = "red") +
      ylab(ylab) +  xlab("Sample size") 
    
  }
  
  if(type== "pBE"){
    if(missing(ylab)) ylab <- "BE Probability"
    se_auc <- sqrt(x$pBE_rep*(1 - x$pBE_rep))
    ymin_auc <- x$pBE_rep - se_auc
    ymax_auc <- x$pBE_rep + se_auc 
    ymax_auc[ymax_auc > 1] <- 1
    ymin_auc[ymin_auc < 0] <- 0
    
    p <- ggplot(x,aes(x=n_subj, y= pBE_rep)) + geom_line(lwd=1.2) +
      geom_errorbar(aes(ymin = ymin_auc,
                        ymax = ymax_auc,
                        width = 0.3), col = "dodgerblue", alpha = 0.5) +
      ylab(ylab) + xlab("Sample size")
  }
  
  print(p)
}


#' Plot quantiles of simulated data
#' 
#' Generates the median and lower/upper quantiles of the simulated data,
#' stratified by drug.
#' 
#' @param x Output of \code{\link{get_data}}
#' @param ci Confidence interval to plot, default is 0.9 (90%) 
#' per FDA guidance.
#' @param ylab Label for y-axis to override defaults.
#' @return A plot.
#' @author Michael Neely
#' @export
plotQ <- function(x, ci = 0.9){
  quants <- x %>% 
    group_by(time,drug) %>%
    summarize(lo = quantile(conc,(1-ci)/2), median = median(conc), hi = quantile(conc,0.5*ci+.5))
  
  p <- ggplot(quants,aes(x=time,y=median,group=drug,col=drug)) + geom_line(lwd=1) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.1, lty = 2, lwd = 0.3) +
    labs(col = "Drug") + ylab("Concentration") + xlab("Time") #+ ggdark::dark_theme_grey()
    
  
  print(p)
}