
#' Read data file for BE simulations
#' 
#' Read a csv file whose format should be subject id, time, conc, period, drug.
#' 
#' @param name Name of the data file.
#' @return A tibble.
#' @author Michael Neely
#' @export

get_data <- function(name){
  return(read_csv(name,col_names = c("id","time","conc","seq","per","drug"),
                  col_types = "fddfff",
                  na = c("", "NA", "."),
                  skip = 1)
         )
}

#' Calculate AUC
#' 
#' Calculate AUC using trapezoidal approximation from tibble read 
#' by \code{\link{get_data}}.
#' 
#' @param x An appropriate tibble.
#' @return A tibble with id and auc.
#' @author Michael Neely
#' @export

get_AUC <- function(x){
  return(Pmetrics::makeAUC(x,conc~time))
}

#' Extract Cmax
#' 
#' Extract Cmax from tibble read by \code{\link{get_data}}.
#' 
#' @param x An appropriate tibble.
#' @return A tibble with id and cmax.
#' @author Michael Neely
#' @export

get_Cmax <- function(x){
  df <- x %>% group_by(id) %>%
    summarize(cmax = max(conc))
 return(df)
}


