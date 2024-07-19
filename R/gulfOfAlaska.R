#' Benthic Foram Abundances from the Gulf of Alaska Reported by Sharon et al. 2021
#' 
#' Loads the foraminifera abundances collected and published by Sharon et al. (2021),
#'  as used in the simulation analyses published in Belanger and Bapst (2023).

#' @details 
#' This data set contains the absolute abundances (number of specimens identified)
#' for 48 species of benthic foraminifera, across 355 samples, taken from Sharon 
#' et al (2021). These samples were collected from the less than 63 μm size fractions
#' taken from Integrated Ocean Drilling Program Expedition 341 Site U141 and the 
#' co-located jumbo piston core, respectively located at 697 meters and 682 meters
#'  water depth in the Gulf of Alaska (Jaeger et al., 2014).

#' @name gulfOfAlaska

#' @rdname gulfOfAlaska

#' @aliases gulfOfAlaska fullDataTable_GOA DCA1_GOA abundData_GOA

#' @docType data

#' @format This dataset is composed of three objects: 
#' 
#' \describe{
#' \item{fullDataTable_GOA}{The full data table containing sample IDs, DCA-1
#' scores and species abundances.}
#' 
#' \item{DCA1_GOA}{A vector of just the DCA-1 scores for each sample.}
#' 
#' \item{abundData_GOA}{The number of specimens identified as a particular
#'  species for each sample.}}
#' 

#' @source 
#' Belanger and Bapst.
#' 
#' Jaeger et al., 2014
#' 
#' Sharon et al. 2021

# @seealso 


#' @keywords datasets

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
#'     "foram_abundances_forSimulations.txt",
#'     header = TRUE)
#' 
#' DCA1_GOA <- fullDataTable_GOA$DCA1
#' abundData_GOA <- fullDataTable_GOA[,-(1:5)]
#' 
#' save(fullDataTable_GOA, DCA1_GOA, abundData_GOA, file = "data/gulfOfAlaska.Rdata")
#' 
#' }
#' 
NULL