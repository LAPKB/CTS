

#extract results
extract <- function(lmer_res){
  #lmer_res is list of n_samp objects, each containing n_trial lmer objects
  n_samp <- length(lmer_res)
  n_trials <- length(lmer_res[[1]])
  
  get_this_lmer <- function(this_lmer){
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
    
    #mse <- 2*(deltaCI/((sqrt(2/n_subj)*qt(.05,2*n_subj-2,lower.tail=FALSE))))^2
    cv_intra <- sqrt(exp(mse)-1)
    
    return(data.frame(de_ratio = de_ratio, ci_lo = ci_lo, ci_up = ci_up, 
                      cv_intra = cv_intra, pBE_rep = pBE_rep, n_subj = n_subj))
  }
  
  all_res <- lmer_res %>% flatten() %>% 
    map(get_this_lmer) %>% 
    map_dfr(~as_tibble(.)) 
  
  sum_all_res <- all_res %>% group_by(n_subj) %>%
    summarize(across(.cols = !matches("pBE_rep"),mean), pBE_rep = sum(pBE_rep)/n_trials)
    
  return(sum_all_res)

}


