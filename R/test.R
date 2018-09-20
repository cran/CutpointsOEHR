#' @title A simulation data to test cutpointsOEHR
#' @description  A dataframe named 'test'contains simulated (t,d,x,x1). The relationship of log relative hazard and x is set to to quandratic, which results in a U-shaped relationship.
#' @docType data
#' @keywords datasets
#' @name test
#' @usage test
#' @format a dataframe contains 200 rows and 4 variables. The 4 varibles are
#' \describe{
#'  \item{t}{simulated times of developing survival outcomes like deathes, relapes, etc.}
#'  \item{d}{censoring indicator, 1 means that survival outcomes are not observed, 0 means survival outcomes are observed. The censoring proportion is set to be 20 percent. }
#'  \item{x}{a continuous variable which has U-shaped relationship with log relative hazard. }
#'  \item{x1}{a continuous variable which has linear relationship with log relative hazard. }
#' }
NULL
