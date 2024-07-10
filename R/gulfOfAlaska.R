#' Foram Abundances from the Gulf of Alaska Reported by Sharon et al. 2021
#' 
#' Loads the abundance data collected and published by Sharon et al. (2021), as used in the simulation analyses published in Belanger and Bapst (2023).
#' 
#' @examples
#' data(gulfOfAlaska)
#' 
#' ####
#' # (This is not to be run)
#' /dontrun{
#' 
#' # Loading the data files used by Belanger & Bapst 2023
#'     # taken from Sharon et al. supplemental
#' fullDataTable_GOA <- read.table(
#'     "~/workspace/simulating_fossil_indices/data/foram_abundances_forSimulations.txt",
#'     header = TRUE)
#' 
#' DCA1_GOA <- fullDataTable_GOA$DCA1
#' abundData_GOA <- fullDataTable_GOA[,-(1:5)]
#' 
#' save(fullDataTable_GOA, DCA1_GOA, abundData_GOA, file = "data/gulfOfAlaska.Rdata")
#' 
#' }
#' 