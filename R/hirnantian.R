#' Katian-Hirnantian Graptolite Assemblages Reported by Sheets et al. (2016)
#' 
#' This is a data set of relative abundances for planktonic graptolite species 
#' from two sections that cross the Katian-Hirnantian boundary in the late Paleozoic.

#' @details 
#' This data set contains the relative abundances 
#' (proportion of specimens identified)
#' for 43 species of planktonic graptolites, 
#' across 34 samples taken from two outcrops: 
#' (a) a section at Vininni Creek in Nevada, USA, and 
#' (b) a second section at Blackstone River in the Yukon, Canada. 
#' Both of these sections cross the Katian-Hirnantian boundary and the
#' relative abundances reported here were provided as supplementary material
#' with Sheets et al. (2016).

#' @name hirnantian

#' @rdname hirnantian

#' @aliases hirnantian graptCommData graptSampleInfo

#' @docType data

#' @format This dataset is composed of two objects: 
#' 
#' \describe{
#' \item{graptCommData}{A data table of sample IDs and
#'  relative species abundances.}
#' 
#' \item{graptSampleInfo}{A data table of additional locality and geologic age
#'  information related to each sample.}}
#' 

#' @source 
#' Sheets, H. David, Charles E. Mitchell, Michael J. Melchin, Jason Loxton, 
#' Petr Štorch, Kristi L. Carlucci, and Andrew D. Hawkins. 
#' "Graptolite community responses to global climate change and 
#' the Late Ordovician mass extinction." Proceedings of the National 
#' Academy of Sciences 113, no. 30 (2016): 8380-8385.

#' @seealso \code{\link{gulfOfAlaska}}

#' @keywords datasets

#' @examples
#' data(hirnantian)
#' 
#' ########################
#' # # (This is not to be run, just showing how data was loaded)
#' #  
#' # # Sheets et al. community abundance data
#' # graptCommData <- read.csv(
#' #    "grapt_abundances_Sheets_et_al_Vinini&Blackstone_01-09-22.csv"
#' #    , row.names = 1, header = TRUE
#' #    )
#' #
#' # # sample specific info
#' # graptSampleInfo <- read.csv(
#' #    "grapt_siteData_SheetsEtAl.csv",
#' #    row.names = 1, header = TRUE, 
#' #    stringsAsFactors = TRUE
#' #    )  
#' #
#' # save(graptCommData, graptSampleInfo, 
#' #    file = "data/hirnantian.Rdata")
#' 
NULL