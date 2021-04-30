


#' @export
#' 
calc_be <- function(x, design = "all", 
                    n = seq(from = 5, to  = 25, by = 5),
                    n_trials = 50, 
                    subj_min = 5, subj_max = 50, subj_step = 5, 
                    seed = 110321
){
  #extract data
  periods <- unique(x$per)
  n_pers <- length(periods) #each period is a replicate for that drug
  drugs <- unique(x$drug)
  n_drugs <- length(drugs)
  n_curves <- x %>% group_by(per, drug) %>%
    summarize(N = n()) #number of profiles (IDs) per condition 
  n_curves <- n_curves$N[1] #assume the same per condition
  
  n_subj <- seq(from = subj_min, to = subj_max, by = subj_step) #step through sample sizes
  n_samp <- length(n_subj)
  set.seed(seed)
  
  
  if(design == "parallel" | design == "crossover" | design == "all"){
    pBE_para <- up_s <- up_ss <- lo_s <- lo_ss <- rep(0, n_samp)
    
    for (i in 1:n_samp) {
      for (b in 1:n_trials) {
        this_list <- list()
        for(j in 1:n_drugs){
          this_list[[j]] <- x %>% filter(per==1 & drug==drugs[j]) %>%
            slice_sample(n = n_subj[i])
        }
        this_df <- this_list %>% map_dfr(~as_tibble(.))
        
        
        res <- t.test(log10(auc)~drug,this_df,paired=F,conf.level = 0.9)
        this_lo <- 10^res$conf.int[1]
        this_up <- 10^res$conf.int[2]
        
        if (this_lo >= 0.8 & this_up <= 1.25) {
          pBE_para[i] <- pBE_para[i] + 1
        }
        lo_s[i] <- lo_s[i] + this_lo
        lo_ss[i] <- lo_ss[i] + lo_s[i]^2
        up_s[i] <- up_s[i] + this_up
        up_ss[i] <- up_ss[i] + up_s[i]^2
        
      }
    }
    
    pBE_para <- pBE_para / n_trials
    lo_s_para <- lo_s/n_trials
    lo_ss_para <- lo_ss/n_trials
    up_s_para <- up_s/n_trials
    up_ss_para <- up_ss/n_trials
    
    return(list(pBE_para = pBE_para, lo_s_para = lo_s_para, up_s_para = up_s_para))
    
  } #end parallel
  
  if(design == "replicate" | design == "all"){
    
    #ctrl1 <- lmeControl(opt='optim', msMaxIter=1000)
    ctrl <- lmerControl(optimizer ="Nelder_Mead")
    
    ci_lcut=.8; ci_ucut=1.25 # lower and upper cutoffs for CI
    
    #function to simulate for replicate designs
    get_lmer <- function(n,par){
      this_samp <- function(i){
        to_keep <- sample(n_curves,n)
        this_df <- x %>% 
          group_by(drug, per) %>%
          slice(to_keep)
        
        if(par=="auc"){
          target_lmer <- lmer(log10(auc) ~  per + drug + (1 | id), data = this_df, control=ctrl) 
        } else {
          target_lmer <- lmer(log10(cmax) ~  per + drug + (1 | id), data = this_df, control=ctrl) 
        }
        return(target_lmer)
      }
      
      res <-  map(1:n_trials, this_samp)
      return(res)
    } 
    #do the simulation
    res_auc <- map(n_subj,get_lmer,"auc")
    res_cmax <- map(n_subj,get_lmer,"cmax")
    return(list(auc = res_auc, cmax = res_cmax))
  
  } # end replication
  
  
} #end function


