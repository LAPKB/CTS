#' Generates bioequivalence statistics 
#'
#' Generates BE statistics for parallel, cross-over and full replication designs.
#' From simulated data with 1000 subjects per drug, the function first calculates 
#' area under the concentration-time curve (AUC) and maximum concentration (Cmax) 
#' for each simulated profiles. Periods are defined as >1 for each replicated
#' administration (e.g., when intra-occasion varibility is simulated.) 
#' Linear regression (\code{\link{lm}}) is used on \code{log10(auc)} and \code{log10(cmax)}
#' for parallel and cross-over designs, and mixed effect modeling (\code{\link{lmer}})
#' is used for the replicated design. 
#'
#' @title Calculation of BE statistics
#' @param data A dataframe or tibble  with these columns:
#' subject id, time, concentrations, period, drug.
#' @param n_trials Integer coding the number of random trials per sample size.
#' @param subj_min Integer defining the minimum number of subjects per trial
#' @param subj_max Integer defining the maximum number of subjects per trial
#' @param subj_step Integer defining the step size between \code{subj_min} and
#' \code{subj_max}, so that the program will test sample sizes 
#' \emph{N = \{subj_min, subj_min + subj_step, subj_min + 2*subj_step,...,subj_max\}}.
#' @param seed Random number seed.
#' @return A list of lists: \[\[1:n_samp\]\]\[\[1:\code{n_trials}\]\]\[\[$par, $cross, 
#' $rep\]\]\[\[1,2]\], where n_samp is the number of sample sizes tested, e.g., 
#' (\code{subj_max} - \code{subj-min}) / \code{subj_step}.  The first item in the 
#' design is auc, and the second is cmax.  For example the comparison of auc using 
#' the third trial of the second sample size for the replicated design would be 
#' x\[\[2\]\]\[\[3\]\]$rep\[\[1\]\].
#' @seealso \code{\link{extractBE}}, \code{\link{plotBE}}
#' @author Michael Neely
#' @export


calc_be <- function(data, 
                    n_trials = 50, 
                    subj_min = 5, subj_max = 50, subj_step = 5, 
                    seed = -17, ci_lcut=.8, ci_ucut=1.25){
  
  #make auc/cmax object from the data
  data2 <- data %>% mutate(id2 = id, id = interaction(id2,per,drug))
  aucV <- get_AUC(data2)$tau
  cmaxV <- get_Cmax(data2)$cmax
  data3 <- data2 %>% group_by(id) %>% distinct(id,.keep_all=T) %>%
    mutate(id=id2) %>% add_column(auc = aucV, cmax = cmaxV) %>%
    select(id,auc,cmax,per,drug)
  
  #extract data
  periods <- unique(data3$per)
  n_pers <- length(periods) #each period is a replicate for that drug
  drugs <- unique(data3$drug)
  n_drugs <- length(drugs)
  n_curves <- data3 %>% group_by(per, drug) %>%
    summarize(N = n()) #number of profiles (IDs) per condition 
  n_curves <- n_curves$N[1] #assume the same per condition
  
  n_subj <- seq(from = subj_min, to = subj_max, by = subj_step) #step through sample sizes
  n_samp <- length(n_subj)
  set.seed(seed)
  
  
  #simulate
  
  ctrl <- lme4::lmerControl(optimizer ="Nelder_Mead")
  
  get_lmer <- function(n){
    #set up progress bar
    pb <- progressr::progressor(along = 1:n_trials, message = paste0("Sample size = ", n ))
    this_samp <- function(i){
      #create subset with n random individuals by drug for parallel design
      par_df <- data3 %>%
        group_by(drug) %>%
        filter(per=="1") %>%
        slice_sample(n=n)
      
      #create subset with the same n random individuals by drug, period
      #for replicate (all periods) and cross-over (period 1 only)
      to_keep <- sample(n_curves,n) #replace=F by default  
      rep_df <- data3 %>% 
        group_by(drug, per) %>%
        slice(to_keep) 
      
      cross_df <- rep_df %>% filter(per=="1" )
      
      lmPar <- list(rep(NA,2)) #for parallel
      lmCross <- list(rep(NA,2)) #for crossover
      lmerRep <- list(rep(NA,2)) #for replicate
      lmPar[[1]] <- suppressWarnings(lm(log10(auc) ~  drug, data = par_df)) 
      lmPar[[2]] <- suppressWarnings(lm(log10(cmax) ~  drug, data = par_df)) 
      lmCross[[1]] <- suppressWarnings(lm(log10(auc) ~  drug, data = cross_df)) 
      lmCross[[2]] <- suppressWarnings(lm(log10(cmax) ~  drug, data = cross_df)) 
      if(n_pers>1) {
        lmerRep[[1]] <- suppressWarnings(lmerTest::lmer(log10(auc) ~  per + drug + (1 | id), data = rep_df, control=ctrl)) 
        lmerRep[[2]] <- suppressWarnings(lmerTest::lmer(log10(cmax) ~  per + drug + (1 | id), data = rep_df, control=ctrl)) 
      } 
      pb()
      
      return(list(par = lmPar, cross = lmCross, rep = lmerRep))
    }
    
    res <- map(1:n_trials, this_samp)
    
    return(res)
    
  } 
  #do the simulation
  res_list <- map(n_subj,get_lmer)
  #[[1:n_samp]][[1:n_trials]][[$par, $cross, $rep]][[auc, cmax]]
  
  #return all objects
  return(res_list)
  
} #end function


