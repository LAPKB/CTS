
#' Read data file for BE simulations
#' 
#' Read a csv file whose format should be subject id, time, conc, period, drug.
#' 
#' @param name Name of the data file.
#' @return A tibble.
#' @author Michael Neely
#' @export

get_data <- function(name){
  return(read_csv(name,col_names = c("id","time","conc","per","drug"),
                  col_types = "fddff",
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

#' Add vertical reference line to plotly
#' 
#' @param x X-intercept
#' @param color Color
#' @param dash Character setting the dash style of lines. Set to a dash type string (solid,
#'  dot, dash, longdash, dashdot, or longdashdot) or a dash length list in px 
#'  (eg 5px,10px,2px,2px).
#' @return A line
#' @author Michael Neely
#' @export

vline <- function(x = 0, color = "red", dash = "solid") {
  list(
    type = "line",
    y0 = 0, 
    y1 = 1, 
    yref = "paper",
    x0 = x, 
    x1 = x, 
    line = list(color = color, dash = dash)
  )
}

#' Add horizontal reference line to plotly
#' 
#' @param y Y-intercept
#' @param color Color
#' @param dash Character setting the dash style of lines. Set to a dash type string (solid,
#'  dot, dash, longdash, dashdot, or longdashdot) or a dash length list in px 
#'  (eg 5px,10px,2px,2px).
#' @return A line
#' @author Michael Neely
#' @export

hline <- function(y = 0, color = "red", dash = "solid") {
  list(
    type = "line",
    x0 = 0, 
    x1 = 1, 
    xref = "paper",
    y0 = y, 
    y1 = y, 
    line = list(color = color, dash = dash)
  )
}