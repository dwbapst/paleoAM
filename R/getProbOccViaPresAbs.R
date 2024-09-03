#' Estimate the Per-Species Probability of Occurrence as a Function of an Environmental Gradient from Presence-Absence Data
#' 
#' This function calculates how species presence (not abundance) varies as a function of an underlying environmental gradient, using only the pattern of presence-absence observed in a given data set (usually data on species abundance). 

#' @details
#' The rationale for this function was that simulating assemblages using the KDEs of species abundance alone (calculated with \code{getSpeciesSpecificRescaledKDE} tended to create very diversity-high assemblages that did not reflect reality. Accounting for species presence at all using separate models brings simulated assemblages much closer to real assemblages.

#' @inheritParams getSpeciesSpecificRescaledKDE

#' @param occurrenceFloor The minimum occurrence for every species, in every bin. 
#' The default is zero -- increasing this value means every species has a 
#' non-zero chance of occurring in every bin, which tends to result in wildly 
#' more diverse assemblages than what is observed in the fossil record, so use with caution.

#' @return
#' An approximate function, created with \code{approx} that describes the 
#' relationship between gradient and probability of occurrence for all species.

# @aliases

#' @seealso
#' \code{\link{getSpeciesSpecificRescaledKDE}}

# @references

#' @examples
#' 
#' # load data
#' data(gulfOfAlaska)
#' 
#' probSpeciesOccur <- getProbOccViaPresAbs(
#'    gradientOrigDCA = DCA1_GOA, 
#'    origAbundData = abundData_GOA
#'    )
#'

#' @name getProbOccViaPresAbs
#' @rdname getProbOccViaPresAbs
#' @export
getProbOccViaPresAbs <- function(
        origAbundData, 
        gradientOrigDCA, 
        occurrenceFloor = 0,
        nBreaksGradientHist = 20
        ){

    # get Probability of Occurrence Along Gradient Using Presence-Absence
    # Get the presence-absence version of the abundance data
    
    # convert abundance data to presence data
    origPresenceData <- origAbundData
    origPresenceData[origPresenceData > 0] <- 1
    # dim(origPresenceData) == dim(origAbundData)
    
    # get histBreaks and histSampleCounts
    gradientHist <- hist(
        gradientOrigDCA, 
        breaks = nBreaksGradientHist, 
        plot = FALSE)
    
    # Count the number of presences across the gradient.
    gradientRepPresence <- apply(
        origPresenceData, 2,
        function(x) gradientOrigDCA[x > 0]                         
        )
    
    # count number of occurrence in each hist bin
    gradientCountPresence <- lapply(gradientRepPresence, function(x) 
        hist(x, breaks = gradientHist$breaks, plot = FALSE)$counts 
        )
    
    # Need to correct for uneven sampling across gradient
        # divide by total number of samples in that bin
    # if occurrenceFloor is not zero (the default)
        # raise floor so every species has a tiny chance of occurring in each bin
    gradientPropPresence <- lapply(gradientCountPresence, function(x)
        (x + occurrenceFloor) / (gradientHist$counts + occurrenceFloor)
        )
    
    # Define the relationship between gradient and probability of occurrence
    # as an approximate function with `approx`
    probSpeciesOccur <- function(gradientValue){
        approxProbs <- lapply(gradientPropPresence, function(x)
            approx(y = x, x = gradientHist$mids, xout = gradientValue)$y
            )  
        }

    # testing lines for checking out calculations
        # when debugging
    #which(probSpeciesOccur(-0.5))
    #probSpeciesOccur(1)
    #probSpeciesOccur(2)
    
    return(probSpeciesOccur)
    }

