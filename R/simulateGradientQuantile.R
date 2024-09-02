#' Repeatedly Simulate Sampled Assemblages at some Gradient Value, and Return a Quantile Based on Recovered Gradient Values
#' 
#' This function simulates assemblages at a single given gradient value, and returns a specified quantile on the recovered gradient values for the sake of defining an envelope around on recovered gradient values.

#' @details
#' This function is most useful with applications like \code{getRecoveredTransitionDuration} which use envelope values to define features of a recovered sequence of gradient values for comparing simulated and empirical data.

#' @inheritParams getTimestepAbundances 
#' @inheritParams sampleFossilSeries
#' @inheritParams quantifyCommunityEcology

#' @param quantileProbs The quantiles to return on the recovered gradient values from the simulated assemblages.

#' @param nSamplesSim The number of samples to simulate.

#' @param gradientValue The gradient value to simulate assemblages at.

#' @return
#' A value for each quantile specified in \code{quantileProbs}. May be multiple values if \code{quantileProbs} is a vector with more than one value.

# @aliases

# @seealso

# @references

# @examples

#' @name simulateGradientQuantile
#' @rdname simulateGradientQuantile
#' @export
simulateGradientQuantile <- function(
            quantileProbs = c(0.95),
            nSamplesSim,
            gradientValue,
            origAbundData,
            kdeRescaled,
            probSpeciesOccur,
            powerRootTransform = 1,
            specimensPerTimestep = 10000,
            nSpecimens
            ){
    
    # previously named as generateBackgroundOnlySimulationQuantile
    
    # generate simulations based only on the background gradient value
        # use these to estimate quantiles which we use to denote 'significance'
        # for the purpose of testing whether an event sample was significant different
        # from the background values
    
    # Figure out which samples will be taken and at what core depths / outcrop heights 
    # First, make a matrix with age, depth, gradient values
    simTimeVar <- data.frame(
        timestep = 1:(nSamplesSim*3),
        # reverse core depth so oldest time is at bottom (biggest depths)
          # 07-08-21: subtract one to start from 0
        sedColumnDepth = ((1:nSamplesSim) -1) * 3,
        gradientValue = gradientValue
        )
        
    # Simulate actual specimen abundances for each time step
    # using KDEs and occurrence probabilities from empirical data
    timestepAbundances <- getTimestepAbundances(
        gradientValues = simTimeVar$gradientValue,
        kdeRescaled = kdeRescaled,
        probSpeciesOccur = probSpeciesOccur,
        specimensPerTimestep = specimensPerTimestep
        )
    colnames(timestepAbundances) <- names(kdeRescaled)
        
    # Running the Simulation
    fossilSeries <- sampleFossilSeries(    
        bioturbIntensity = 0, 
        bioturbZoneDepth = 0,  
        distBetweenSamples = 0, 
        sampleWidth = 3, 
        simTimeVar = simTimeVar,
        timestepAbundances = timestepAbundances,
        nSpecimens = nSpecimens
        )
        
    # get sample DCA1 scores
    ecologyOutList <- quantifyCommunityEcology(
        origAbundData = origAbundData,
        fossilSeries = fossilSeries,
        singularDCA = TRUE,
        inclusiveDCA = FALSE,
        rawDCA = FALSE,
        powerRootTransform = powerRootTransform 
        )

    DCA1 <- ecologyOutList$scoreDCA1_singular
    
    quantileOut <- quantile(
        x = DCA1, probs = quantileProbs, na.rm = TRUE
        )
    
    return(quantileOut)
    }
