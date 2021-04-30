#' Clinical trial simlation for BE studies
#' 
#' This package contains functions to generate statistics for various sample
#' sizes and trial designs based on models created by PK-SIM + Pmetrics.
#'
#' @name CTS-package
#' @aliases CTS
#' @docType package
#' @author Michael Neely MD
#' \url{http://www.lapk.org}
#' @keywords package
#'
#' @importFrom readr read_csv
#' @importFrom Pmetrics makeAUC
#' @importFrom nlme lmeControl lme
#' @importFrom ggplot2 ggplot
#' #' @importFrom dplyr select arrange filter mutate transmute group_by row_number
#' @importFrom foreach %dopar%
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_polygon 
#' scale_x_log10 scale_x_continuous scale_y_log10 scale_y_continuous xlab ylab
#' theme ggtitle element_blank
#' @importFrom purrr map reduce map_chr keep pluck %>% 
#' @importFrom tibble as_tibble
#' @importFrom tidyr pivot_longer pivot_wider unnest extract separate fill
#' @importFrom grDevices col2rgb dev.off devAskNewPage gray.colors jpeg 
#' pdf png postscript rgb setEPS
#' @importFrom graphics abline arrows axTicks axis boxplot hist legend lines 
#' par plot points polygon rect rug segments text
#' @importFrom stats aggregate anova approx as.formula binom.test coef 
#' complete.cases confint cor cor.test cov cov.wt cov2cor density dnorm 
#' get_all_vars glm kmeans kruskal.test ks.test
#' lm median model.frame pchisq pnorm predict
#' pt qchisq qnorm qqline qqnorm qqplot qt
#' quantile rnorm runif sd shapiro.test step
#' t.test terms time var weighted.mean wilcox.test
#' @importFrom utils compareVersion data flush.console glob2rx head 
#' install.packages news packageVersion read.table setTxtProgressBar str
#' tail txtProgressBar write.csv write.table

NULL

