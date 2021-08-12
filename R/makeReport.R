#' Generate BE report
#' 
#' @param sim Raw simulated data from \code{\link{get_data}}.
#' @param be Output of \code{\link{extractBE}}
#' @return HTML document will open in browser
#' @author Michael Neely
#' @export

make_report <- function(sim, be){
  
  outputFile <- "BEreport.html"
  reportFile <- system.file("report", "BEreport.Rmd", package = "CTS")
  rmarkdown::render(reportFile, params = list(sim=sim, be=be), 
                    output_file = outputFile, output_dir = getwd())
  browseURL(outputFile)
  return(invisible())
}