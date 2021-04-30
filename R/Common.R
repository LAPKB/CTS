#tidyverse


#read data
#data format should be subject, time, conc, sequence, period, drug
#expected number 1000 per condition, with 2 sequences, and 4 periods
#total is 8000 concentrations per time point

#' @export
get_data <- function(name){
  return(read_csv(name,col_names = c("id","time","conc","seq","per","drug"),
                  col_types = "fddfff",
                  na = c("", "NA", "."),
                  skip = 1)
         )
}

#calculate AUC

#' @export
get_auc <- function(x){
  return(makeAUC(x,conc~time))
}

#get Cmax

#' @export
get_cmax <- function(x){
 return(tapply(x$conc,x$id,max))

}


#Test BE for parallel design

# '@export
test_parallel <- function(stat,drug,paired){
  t_res <- t.test(log10(stat)~drug, paired = paired, conf.level = 0.9)
  return(10**t_res$conf.int)
}
