#' @details

#' @inheritParams

#' @param

#' @return

#' @aliases

#' @seealso

#' @references

#' @examples

#' @param origAbundData 
#'
#' @param gradientOrigDCA 
#' @param occurrenceFloor 
#' @param nBreaksGradientHist 
#' @param plot 
#'
#' @name
#' @rdname
#' @export

getProbOccViaPresAbs <- function(
        origAbundData, 
        gradientOrigDCA, 
        occurrenceFloor = 0,
        nBreaksGradientHist = 20, 
        plot = FALSE
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
    
    if(plot){
        hist(gradientRepPresence$EpistominellaPacifica, 
             xlab = "Gradient (DCA 1 Score)", 
            breaks = gradientHist$breaks,
             main = "Samples where Epistominella pacifica occurs")
        hist(gradientRepPresence$LagenaSpp, 
             xlab = "Gradient (DCA 1 Score)", 
            breaks = gradientHist$breaks,
             main = "Samples where Lagena Spp occurs")
        }
    
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
    
    if(plot){
        plot(gradientHist$mids, 
            gradientPropPresence$EpistominellaPacifica, 
             xlab = "Gradient (DCA 1 Score)",
            type = "b",
             main = "Proportion of Samples with Epistominella pacifica sampled")
        plot(gradientHist$mids, 
            gradientPropPresence$LagenaSpp, 
             xlab = "Gradient (DCA 1 Score)", 
            type = "b",
             main = "Proportion of Samples with Lagena Spp sampled")
        }
    
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

