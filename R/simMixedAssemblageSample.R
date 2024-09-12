#' Simulate a Mixed Fossil Assemblage Composed of Communities at Different Gradient Values, and Sample the Lumped Assemblage
#' 
#' This function simulate a mixed fossil assemblage by simulating a series of communities across a defined range of gradient values, lumps them into a single mixed assemblage, and then samples that assemblage as defined by the user.

#' @details
#' This function is mainly written for simulating what artificial mixtures of assemblages at different gradient values would look like if sampled and assumed to be a single cohesive assemblage.

#' @inheritParams getTimestepAbundances 

#' @inheritParams sampleFossilSeries

#' @return
#' A matrix containing the species abundances in the resulting mixed assemblage.

# @aliases

# @seealso

# @examples


#' @name simMixedAssemblageSample
#' @rdname simMixedAssemblageSample
#' @export
simMixedAssemblageSample <- function(
            kdeRescaled, 
            probSpeciesOccur, 
            gradientValues, 
            specimensPerTimestep,
            nSpecimens     
            ){
    # simulate a mixed sample, generated from lumping 'background'
        # and event assemblages into a single assemblage    

    # lump, then pick the sample
    lumpedSample <- getTimestepAbundances(
            kdeRescaled = kdeRescaled, 
            probSpeciesOccur = probSpeciesOccur, 
            gradientValues = gradientValues, 
            specimensPerTimestep = specimensPerTimestep
            )

    lumpedSample <- colSums(lumpedSample)    
    
    # number of species = length(kdeRescaled)
    
    # un-table() the lumped community abundance data
    lumpedSample <- rep(1:length(kdeRescaled), lumpedSample)
    
    # down-sample the lumped specimens to "nSpecimens" 
    pickedSample <- sample(x = lumpedSample, 
               replace = FALSE, size = nSpecimens)
    
    # down-sample the lumped specimens to "nSpecimens" 
    pickedSample <- tabulate(pickedSample, nbins = length(kdeRescaled))
        
    return(pickedSample)
    }
