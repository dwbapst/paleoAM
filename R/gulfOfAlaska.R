#' Foram Abundances from the Gulf of Alaska Reported by Sharon et al. 2021
#' 
#' @examples
#' data(gulfOfAlaska)
#' 
#' ####
#' # (This is not to be run)
#' /dontrun{
#' 
#' # Loading the data files from Belanger & Bapst suppmat
#' 
#' 
#' 
#' 
#' 
#' }
#' 

fullDataTable_GOA <- read.table(
    "~/workspace/simulating_fossil_indices/data/foram_abundances_forSimulations.txt",
    header = TRUE)

DCA1_GOA <- fullDataTable_GOA$DCA1
abundData_GOA <- fullDataTable_GOA[,-(1:5)]

save(fullDataTable_GOA, DCA1_GOA, abundData_GOA, file = "data/gulfOfAlaska.Rdata")
