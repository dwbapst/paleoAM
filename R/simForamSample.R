simForamSample <- function(
            kdeRescaled, 
            sampleSpeciesGradient, 
            gradientValues, 
            specimensPerTimestep,
            nSpecies,
            nPickedSpecimens     
            ){
    # simulate a mixed sample, generated from lumping 'background'
        # and event assemblages into a single assemblage    

    # lump, then pick the sample
    lumpedSample <- colSums(
        getTimestepAbundances(
            kdeRescaled = kdeRescaled, 
            sampleSpeciesGradient = sampleSpeciesGradient, 
            gradientValues = gradientValues, 
            specimensPerTimestep = specimensPerTimestep
            )
        )
    
    # un-table() the lumped community abundance data
    lumpedSample <- rep(1:nSpecies, lumpedSample)
    
    # down-sample the lumped specimens to "nSpecimensPicked" 
    pickedSample <- sample(x = lumpedSample, 
               replace = FALSE, size = nPickedSpecimens)
    
    # down-sample the lumped specimens to "nSpecimensPicked" 
    pickedSample <- tabulate(pickedSample, nbins = nSpecies)    
        
    return(pickedSample)
    }
