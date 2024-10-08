#' Sample Fossil Assemblage Series
#'
#' Given a time-series of 'true' fossil assemblages simulated in precise time,
#' this function then chunks that 'true' ecological signal into sedimentary
#' packages, which contain specimens from assemblages spanning the time interval
#' during which that sediment accumulated. Further more, the inclusion of 
#' specimens from even more distant assemblages is used to model bioturbation.

#' @details
#' This function is where bioturbation processes are handled, 
#' as well as time-averaging from samples capturing several 
#' sedimentary horizons reflecting multiple original fossil assemblages.
#' 
#' This function is generally run after running \code{\link{getTimestepAbundances}}. 
#' Most users will likely never run either function, instead running
#' \code{\link{simulateFossilAssemblageSeries}}.


#' @inheritParams calculateImplicitParameters

#' @param bioturbIntensity The degree of mixing within the bioturbation zone, as a value between 0 and 1. When intensity is 1, a given sample will consist only 
 
#' @param bioturbZoneDepth The sediment depth to which bioturbation occurs. For example, Bioturbation depth varies considerably in the modern ocean, but is often around 10 centimeters -- with the top ten centimeters of sediment (and the organic remains in those ten centimeters of sediment) being regularly moved up and down by organism activity. For the purposes of this model, a bioturbation zone depth of 10 centimeters means that sampling a centimeter of sediment at location X, the apparent fossil assemblage that would be recovered is just as likely to include specimens that were deposited five centimeters away as those deposited at location X.

#' @param nSpecimens The number of specimens selected in each individual sample.

#' @param distBetweenSamples The sedimentary thickness between successive 
#' samples, in the same units as \code{sampleWidth}.
 
#' @param simTimeVar A data-frame specifying time-steps, sedimentary depth and environmental gradient values for simulating a time-series of sampled fossil assemblages.

#' @param timestepAbundances A matrix containing abundances for species as a series of simulated assemblages, output by \code{\link{getTimestepAbundances}}.
 


#' @return
#' A list composed of four components:
#' \code{simTimeVar}, the input data-frame specifying time-steps, sedimentary depth and environmental gradient values;
#' \code{abundanceTable}, a table of the abundances of species in each sample;
#' \code{sampleIntervals}, a table specifying when in time each sample 'begins' and 'ends' in time (based on the sedimentation rate),
#' and \code{bioturbIntervals}, a table specifying which intervals are 'included' in a sample 

# @aliases

#' @seealso
#' This function is generally run after running
#' \code{\link{getTimestepAbundances}}. Most users will likely never run either function, instead running \code{\link{simulateFossilAssemblageSeries}}.

# @references

# @examples


#' @name sampleFossilSeries
#' @rdname sampleFossilSeries
#' @export
sampleFossilSeries <- function(
            bioturbIntensity, 
            bioturbZoneDepth,  
            distBetweenSamples, 
            sampleWidth, 
            simTimeVar, 
            timestepAbundances,
            nSpecimens
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
    while(max(sampleIntervals[,1]) < (
            max(simTimeVar$sedColumnDepth) - sampleWidth - distBetweenSamples
            )){
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
        bioturbIntervals[bioturbIntervals[,1] < min(simTimeVar$sedColumnDepth), 1
            ] <- min(simTimeVar$sedColumnDepth)
        bioturbIntervals[bioturbIntervals[,2] < min(simTimeVar$sedColumnDepth), 2
            ] <- min(simTimeVar$sedColumnDepth)
        bioturbIntervals[bioturbIntervals[,1] > max(simTimeVar$sedColumnDepth), 1
            ] <- max(simTimeVar$sedColumnDepth)
        bioturbIntervals[bioturbIntervals[,2] > max(simTimeVar$sedColumnDepth), 2 
            ] <- max(simTimeVar$sedColumnDepth)
        
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
            timesteps_After   <- (simTimeVar$sedColumnDepth <= bioturbIntervals[i,1])
            timesteps_Before  <- (bioturbIntervals[i,2] <= simTimeVar$sedColumnDepth)
            selectedTimesteps <- timesteps_After & timesteps_Before
            selectedTimesteps <- which(selectedTimesteps)

            # 'sample timesteps' are those 'selected timesteps' 
                # that are the timesteps in the 'actual' sample
            timesteps_After  <- (
                simTimeVar$sedColumnDepth[selectedTimesteps] <= sampleIntervals[i,1])
            timesteps_Before <- (
                sampleIntervals[i,2] <= simTimeVar$sedColumnDepth[selectedTimesteps])
            
            sampleTimesteps  <- timesteps_After & timesteps_Before
            sampleTimesteps  <- which(sampleTimesteps)
            
            nSampleSpec <- round((1 - bioturbProportion) * nSpecimens)   
            nBioturbSpec <- round(bioturbProportion * nSpecimens)
            
        }else{
            # which timestep-layers are in this sample
            timesteps_After   <- (simTimeVar$sedColumnDepth <= sampleIntervals[i,1])
            timesteps_Before  <- (sampleIntervals[i,2] <= simTimeVar$sedColumnDepth)
            selectedTimesteps <- timesteps_After & timesteps_Before
            selectedTimesteps <- which(selectedTimesteps)
            
            sampleTimesteps <- 1:length(selectedTimesteps)
            
            nSampleSpec <- nSpecimens
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
        # down-sample the lumped specimens to "nSpecimens" 
        pickedSample <- sample(x = lumpedSample, 
               replace = FALSE, size = nSampleSpec)
                
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
                                    replace = FALSE, size = nBioturbSpec)
            # combine
            pickedSample <- c(pickedSample, pickedBioturb)
            }
        
        # down-sample the lumped specimens to "nSpecimens" 
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
