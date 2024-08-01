#' @details

#' @inheritParams

#' @param

#' @return

#' @aliases

#' @seealso

#' @references

#' @examples

#' @param kdeRescaled 

#' @param probSpeciesOccur 

#' @param gradientValues 

#' @param specimensPerTimestep 

#' @param nSpecies 

#' @param nPickedSpecimens 

#' @name
#' @rdname
#' @export


simAssemblageSample <- function(
            kdeRescaled, 
            probSpeciesOccur, 
            gradientValues, 
            specimensPerTimestep,
            nSpecies,
            nPickedSpecimens     
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
        
    # un-table() the lumped community abundance data
    lumpedSample <- rep(1:nSpecies, lumpedSample)
    
    # down-sample the lumped specimens to "nSpecimensPicked" 
    pickedSample <- sample(x = lumpedSample, 
               replace = FALSE, size = nPickedSpecimens)
    
    # down-sample the lumped specimens to "nSpecimensPicked" 
    pickedSample <- tabulate(pickedSample, nbins = nSpecies)    
        
    return(pickedSample)
    }
