#' Benthic Foram Abundances from the Gulf of Alaska Reported by Sharon et al 2021
#' 
#' Loads the foraminifera abundances collected and published by Sharon et al. (2021),
#'  as used in the simulation analyses published in Belanger and Bapst (2023).

#' @details 
#' This data set contains the absolute abundances (number of specimens identified)
#' for 48 species of benthic foraminifera, across 355 samples, taken from Sharon 
#' et al (2021). These samples were collected from the less than 63 micrometer size-fractions
#' taken from Integrated Ocean Drilling Program Expedition 341 Site U141 and the 
#' co-located jumbo piston core, respectively located at 697 meters and 682 meters
#' water depth in the Gulf of Alaska (Jaeger et al., 2014).

#' @name gulfOfAlaska

#' @rdname gulfOfAlaska

#' @aliases gulfOfAlaska fullDataTable_GOA DCA1_GOA abundData_GOA

#' @docType data

#' @format This data set is composed of three objects: 
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
#' Belanger, Christina L., and David W. Bapst. 2023.
#' "Simulating our ability to accurately detect abrupt 
#' changes in assemblage-based paleoenvironmental proxies." 
#' Palaeontologia Electronica 26 (2), 1-32
#' 
#' Jaeger, J.M., Gulick, S.P.S., LeVay, L.J., Asahi, H., Bahlburg, H., 
#' Belanger, C.L., Berbel, G.B.B., Childress, L.B., Cowan, E.A., Drab, L., 
#' Forwick, M., Fukumura, A., Ge, S., Gupta, S.M., Kioka, A., Konno, S., 
#' März, C.E., Matsuzaki, K.M., McClymont, E.L., Mix, A.C., Moy, C.M., 
#' Müller, J., Nakamura, A., Ojima, T., Ridgway, K.D., Rodrigues Ribeiro, F., 
#' Romero, O.E., Slagle, A.L., Stoner, J.S., St-Onge, G., Suto, I., 
#' Walczak, M.H., and Worthington, L.L., 2014. Expedition 341 summary. 
#' In Jaeger, J.M., Gulick, S.P.S., LeVay, L.J., and the 
#' Expedition 341 Scientists, Proceedings of IODP, 341: College Station, TX 
#' (Integrated Ocean Drilling Program). 
#' 
#' Sharon, Christina Belanger, Jianghui Du, and Alan Mix. 
#' "Reconstructing paleo‐oxygenation for the last 54,000 years in the 
#' Gulf of Alaska using cross‐validated benthic foraminiferal and 
#' geochemical records." Paleoceanography and Paleoclimatology 
#' 36, no. 2 (2021): e2020PA003986.

#' @seealso \code{\link{hirnantian}}


#' @keywords datasets

#' @examples
#' data(gulfOfAlaska)
#' 
#' ########################
#' # (This is not to be run, just showing how data was loaded)
#' #
#' # # Loading the data files used by Belanger & Bapst 2023
#' #     # taken from Sharon et al. supplemental
#' # fullDataTable_GOA <- read.table(
#' #    "foram_abundances_forSimulations.txt",
#' #    header = TRUE)
#' # 
#' # DCA1_GOA <- fullDataTable_GOA$DCA1
#' # abundData_GOA <- fullDataTable_GOA[,-(1:5)]
#' #
#' # save(fullDataTable_GOA, DCA1_GOA, abundData_GOA, 
#' #    file = "data/gulfOfAlaska.Rdata")
#' 
#' 
NULL