#' Extraxt bioequivalence statistics 
#'
#' Extracts BE statistics from results of (\code{\link{calc_BE}}). 
#'
#' @title Extract  BE statistics
#' @param lmer_res The output of (\code{\link{calc_BE}}).
#' @param ci Confidence interval for geometric mean ratio, default is 0.9 (90%) 
#' per FDA guidance.
#' @param ci_lcut Lower confidence interval cutoff for geometric mean AUC or Cmax to be
#' considered BE, default 0.8 per FDA guidance.
#' @param ci_ucut Upper confidence interval cutoff for geometric mean AUC or Cmax to be
#' considered BE, default 1.25 per FDA guidance.
#' @return A list of lists: \[\[1:n_samp\]\]\[\[1:\code{n_trials}\]\]\[\[$par, $cross, 
#' $rep\]\]\[\[1,2]\], where n_samp is the number of sample sizes tested, e.g., 
#' (\code{subj_max} - \code{subj-min}) / \code{subj_step}.  The first item in the 
#' design is auc, and the second is cmax.  For example the comparison of auc using 
#' the third trial of the second sample size for the replicated design would be 
#' x\[\[2\]\]\[\[3\]\]$rep\[\[1\]\].
#' @seealso \code{\link{extractBE}}, \code{\link{plotBE}}
#' @author Michael Neely
#' @export
extract_be <- function(lmer_res, ci = 0.9, ci_lcut = 0.8, ci_ucut = 1.25){
  
  n_samp <- length(lmer_res)
  n_trials <- length(lmer_res[[1]])
  rep_design <- inherits(lmer_res[[1]][[1]]$rep[[1]],"lmerModLmerTest")
  
  get_summary <- function(design,index){
    pb <- progressr::progressor(along = 1:length(flatten(lmer_res)), message = 
                       paste0("Extracting ",c("Parallel ", "Cross-over ","Replicate ")[design], 
                              c("AUC","Cmax")[index]))
    
    get_this_lmer <- function(this_lmer,design){
      if(design <= 2){ #parallel or cross over
        sum_lm <- summary(this_lmer)
        n_subj <- length(this_lmer$residuals)/length(this_lmer$xlevels$drug)
        dd_eff <- sum_lm$coefficients[2,1] #drug effect
        se <- sum_lm$coefficients[2,2] #standard error
        df <- sum_lm$df[2] #degrees of freedom
        de_ratio <- exp(dd_eff)
        deltaCI <- qt((1-ci)/2,df,lower.tail=FALSE)*se
        ci_lo <- de_ratio*exp(-deltaCI)
        ci_up <- de_ratio*exp(deltaCI)
        if ((ci_lo>=ci_lcut) & (ci_up<=ci_ucut)) {pBE_rep <- 1} else {pBE_rep <- 0}
        pb()
        
        return(data.frame(de_ratio = de_ratio, ci_lo = ci_lo, ci_up = ci_up, 
                          pBE_rep = pBE_rep, n_subj = n_subj))
      } 
      
      if(design == 3){ #replicate
        # if(!is.na(this_lmer)){
          sum_lmer <- summary(this_lmer)$coefficients
          n_subj <- summary(this_lmer)$ngrps
          dd_eff <- sum_lmer[3,1] #drug effect
          se <- sum_lmer[3,2] #standard error
          df <- sum_lmer[3,3] #degrees of freedom
          de_ratio <- exp(dd_eff)
          deltaCI <- qt(.05,df,lower.tail=FALSE)*se
          ci_lo <- de_ratio*exp(-deltaCI)
          ci_up <- de_ratio*exp(deltaCI)
          if ((ci_lo>=ci_lcut) & (ci_up<=ci_ucut)) {pBE_rep <- 1} else {pBE_rep <- 0}
          
          mse <- 2*(deltaCI/((sqrt(2/n_subj) * qt(.05, 2*n_subj-2, lower.tail=FALSE))))^2
          cv_intra <- sqrt(exp(mse)-1)
        # } else { #replicate is missing
        #   de_ratio <- NA
        #   ci_lo <- NA
        #   ci_up <- NA
        #   cv_intra <- NA
        #   pBE_rep <- NA
        #   n_subj <- NA
        # }
        
        pb()
        
        return(data.frame(de_ratio = de_ratio, ci_lo = ci_lo, ci_up = ci_up, 
                          cv_intra = cv_intra, pBE_rep = pBE_rep, n_subj = n_subj))
      }
      
    }
    
    this_res <- lmer_res %>% flatten() %>%
      map(c(design,index)) %>% 
      #design 1=parallel; 2=cross-over; 3=replicate
      #index 1 = auc, 2 = cmax 
      map(get_this_lmer,design) %>% 
      map_dfr(~as_tibble(.)) %>% 
      group_by(n_subj) %>%
      summarize(across(.cols = !matches("pBE_rep"),mean), pBE_rep = sum(pBE_rep)/n_trials)
  }
  
  auc_par <- get_summary(1,1)
  cmax_par <- get_summary(1,2)
  auc_cross <- get_summary(2,1)
  cmax_cross <- get_summary(2,2)
  if(rep_design){
    auc_rep <- get_summary(3,1)
    cmax_rep <- get_summary(3,2)
  } else {
    auc_rep <- NA
    cmax_rep <- NA
  }

    
  return(list(auc_par = auc_par, cmax_par = cmax_par, 
              auc_cross = auc_cross, cmax_cross = cmax_cross, 
              auc_rep = auc_rep, cmax_rep = cmax_rep))

}


