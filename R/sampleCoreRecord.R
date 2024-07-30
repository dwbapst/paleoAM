# sampleFossilAssemblageSeries

# importantly, this function is where bioturbation is handled

sampleFossilAssemblageSeries <- function(
            bioturbIntensity, 
            bioturbZoneDepth,  
            distBetweenSamples, 
            sampleWidth, 
            simTimeVar, 
            timestepAbundances,
            nSpecimensPicked
            ){
            
    nSpecies <- ncol(timestepAbundances)
    speciesNames <- colnames(timestepAbundances)
    
    # Now work out what samples will be taken, 
        # and what assemblages have to be simulated for each sample
        # based on sampleWidth and bioturbZoneDepth
    #
    # make a matrix
    sampleIntervals <- matrix(c(distBetweenSamples, distBetweenSamples + sampleWidth), 1, 2)
    # run a while loop
    while(max(sampleIntervals[,1]) < (max(simTimeVar$coreDepth) - sampleWidth - distBetweenSamples)){
        newSampleInterval <- sampleIntervals[nrow(sampleIntervals),2] + distBetweenSamples
        newSampleInterval <- c(newSampleInterval, newSampleInterval + sampleWidth)
        sampleIntervals <- rbind(sampleIntervals, newSampleInterval)
        }
    
    #drop all but the last row)
    sampleIntervals <- sampleIntervals[1:(nrow(sampleIntervals) - 1), , drop = FALSE]
    rownames(sampleIntervals) <- NULL
    sampleIntervals <- sampleIntervals[ , 2:1, drop = FALSE]
    colnames(sampleIntervals) <- c("sampleDepth_start", "sampleDepth_end")
    # number of samples
    nSamples <- nrow(sampleIntervals)
    # record the mid-depth for each sample
    sampleMidDepth <- apply(sampleIntervals, 1, mean)
    
    # Get symmetric bioturbation intervals 
    # by working backward from `sampleMidDepth and bioturbZoneDepth
    if((bioturbZoneDepth > sampleWidth) & (bioturbIntensity > 0)){
        
        # what is the proportion of picked specimens that come from the 
        # 'exterior' bioturbated portion?
        bioturbProportion <- bioturbIntensity * (1 - (sampleWidth / bioturbZoneDepth))
        
        # bioturbation going UP (above the sample) and going DOWN (below the sample)
        bioturbIntervals <- sampleMidDepth + bioturbZoneDepth/2
        bioturbIntervals <- cbind(bioturbIntervals, sampleMidDepth - bioturbZoneDepth/2)
        #
        bioturbIntervals[bioturbIntervals[,1] < min(simTimeVar$coreDepth), 1] <- min(simTimeVar$coreDepth)
        bioturbIntervals[bioturbIntervals[,2] < min(simTimeVar$coreDepth), 2] <- min(simTimeVar$coreDepth)
        bioturbIntervals[bioturbIntervals[,1] > max(simTimeVar$coreDepth), 1] <- max(simTimeVar$coreDepth)
        bioturbIntervals[bioturbIntervals[,2] > max(simTimeVar$coreDepth), 2 ] <- max(simTimeVar$coreDepth)
        
    }else{

        # no bioturbation in this simulation
        bioturbProportion <- 0
        bioturbIntervals <- matrix(NA, nrow(sampleIntervals), ncol(sampleIntervals))
        
        }
    
    ## Build Abundance Table for Simulated Picked-Samples
        # Build picked-sample abundance table and fill with species abundances 
        # for each picked sample, looking at each sample individually. 
    # Treat each sample as a lumped mega-sample, 
        # where each timestep contained within a sample is a separately simulated community. 
    # We lump all the sampled specimens from all the timesteps within a picked sample together, 
        # and resample to approximate the time-averaging of the real microfossil record.
        
    # make an empty species by sites matrix for abundances
    abundanceTable <- matrix(0, nrow = nSamples, ncol = nSpecies)
    colnames(abundanceTable) <- speciesNames
    
    # Now run the simulation!  
        
    for(i in 1:nSamples){
        # first determine if bioturbation needs to be considered, 
        # if so, indicate the bioturb intervals 
        if(bioturbProportion > 0){        
            #IF bioturbZoneDepth > sampleWidth
            
            # selected timesteps are the sample timesteps + the timesteps mixed in by bioturbation
            timesteps_After   <- (simTimeVar$coreDepth <= bioturbIntervals[i,1])
            timesteps_Before  <- (bioturbIntervals[i,2] <= simTimeVar$coreDepth)
            selectedTimesteps <- timesteps_After & timesteps_Before
            selectedTimesteps <- which(selectedTimesteps)

            # 'sample timesteps' are those 'selected timesteps' that are the timesteps in the 'actual' sample
            timesteps_After  <- (simTimeVar$coreDepth[selectedTimesteps] <= sampleIntervals[i,1])
            timesteps_Before <- (sampleIntervals[i,2] <= simTimeVar$coreDepth[selectedTimesteps])
            sampleTimesteps  <- timesteps_After & timesteps_Before
            sampleTimesteps  <- which(sampleTimesteps)
            
            nSampleSpecPicked <- round((1 - bioturbProportion) * nSpecimensPicked)   
            nBioturbSpecPicked <- round(bioturbProportion * nSpecimensPicked)
            
        }else{
            # which timestep-layers are in this sample
            timesteps_After   <- (simTimeVar$coreDepth <= sampleIntervals[i,1])
            timesteps_Before  <- (sampleIntervals[i,2] <= simTimeVar$coreDepth)
            selectedTimesteps <- timesteps_After & timesteps_Before
            selectedTimesteps <- which(selectedTimesteps)
            
            sampleTimesteps <- 1:length(selectedTimesteps)
            
            nSampleSpecPicked <- nSpecimensPicked
            }
        
        if(length(selectedTimesteps)<1){
            stop("No sedimentary samples within this interval?!")
            }
        if(length(sampleTimesteps)<1){
            stop("No sedimentary samples within this interval?!")
            }        
        
        # lump the abundances together
        # make an empty matrix
        lumpedAbundanceTable <- timestepAbundances[selectedTimesteps, , drop = FALSE]
        # now do colSums to lump the abundances together for *just* the sample
        lumpedSample <- colSums(lumpedAbundanceTable[sampleTimesteps, , drop = FALSE])
        # un-table() the lumped community abundance data
        lumpedSample <- rep(1:nSpecies, lumpedSample)
        # down-sample the lumped specimens to "nSpecimensPicked" 
        pickedSample <- sample(x = lumpedSample, 
               replace = FALSE, size = nSampleSpecPicked)
                
        # combine with the bioturbated interval *beyond* the sample, if there is bioturbation
        
        # if there are selected timesteps outside of sample timesteps
        if(length(selectedTimesteps) > length(sampleTimesteps)){
            lumpedBioturb <- colSums(lumpedAbundanceTable[-sampleTimesteps, , drop = FALSE])
            # un-table() the community abundance data
            lumpedBioturb <- rep(1:nSpecies, lumpedBioturb)
            
            if(length(lumpedBioturb)<1){
                stop("No species found in surrounding bioturbated layers?")
                }
            
            # and now sample the bioturbated remainder interval
            pickedBioturb <- sample(x = lumpedBioturb, 
                                    replace = FALSE, size = nBioturbSpecPicked)
            # combine
            pickedSample <- c(pickedSample, pickedBioturb)
            }
        
        # down-sample the lumped specimens to "nSpecimensPicked" 
        pickedSample <- tabulate(pickedSample, nbins = nSpecies)
        
        # record in abundance
        abundanceTable[i,] <- pickedSample
        }
    
    fossilSeries <- list(
        simTimeVar = simTimeVar,
        abundanceTable = abundanceTable,
        sampleIntervals = sampleIntervals,
        bioturbIntervals = bioturbIntervals
        )
    
    return(fossilSeries)
    }
