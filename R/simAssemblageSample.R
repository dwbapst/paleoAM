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

#' @param nSpecimens 

#' @name
#' @rdname
#' @export


simAssemblageSample <- function(
            kdeRescaled, 
            probSpeciesOccur, 
            gradientValues, 
            specimensPerTimestep,
            nSpecies,
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
        
    # un-table() the lumped community abundance data
    lumpedSample <- rep(1:nSpecies, lumpedSample)
    
    # down-sample the lumped specimens to "nSpecimens" 
    pickedSample <- sample(x = lumpedSample, 
               replace = FALSE, size = nSpecimens)
    
    # down-sample the lumped specimens to "nSpecimens" 
    pickedSample <- tabulate(pickedSample, nbins = nSpecies)    
        
    return(pickedSample)
    }
