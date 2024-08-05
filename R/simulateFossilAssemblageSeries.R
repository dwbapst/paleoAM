#' @details

#' @inheritParams

#' @param

#' @return

#' @aliases

#' @seealso

#' @references

#' @examples

#' @param kdeRescaled 
#'
#' @param probSpeciesOccur 
#' @param origAbundData 
#' @param eventChangeScale 
#' @param bgGradientValue 
#' @param fullGradientRange 
#' @param eventSampleWidthRatio 
#' @param sampleWidth 
#' @param eventDuration 
#' @param sedRatePerTimestep 
#' @param samplingCompleteness 
#' @param transitionDurationRatio 
#' @param bioturbDepthRatio 
#' @param bioturbIntensity 
#' @param nEvents 
#' @param nSpecimens 
#' @param specimensPerTimestep 
#' @param halfGradientOnly 
#' @param useTransformedRelAbundance 
#' @param projectIntoOrigDCA 
#' @param powerRootTransform 

#' @param singularDCA 
#' @param inclusiveDCA 
#' @param rawDCA 
#' @param includeInitialBackgroundPhase 
#' @param plot 
#'
#' @name
#' @rdname
#' @export


simulateFossilAssemblageSeries <- function(
                kdeRescaled,
                probSpeciesOccur,
                origAbundData,
                eventChangeScale,
                bgGradientValue,
                fullGradientRange,
                eventSampleWidthRatio,
                sampleWidth,
                eventDuration,
                sedRatePerTimestep,
                samplingCompleteness,
                transitionDurationRatio,
                bioturbDepthRatio,
                bioturbIntensity,     
                nEvents,
                nSpecimens,
                #
                specimensPerTimestep = 10000,
                halfGradientOnly = "full",
                useTransformedRelAbundance = TRUE,
                projectIntoOrigDCA = TRUE,
                powerRootTransform = 1,  
                maxSampleTimeStep = 500,
                minSampleTimeStep = 3,
                singularDCA = TRUE,
                inclusiveDCA = FALSE, # BAD IDEA
                rawDCA = FALSE,
                includeInitialBackgroundPhase = FALSE,
                # runChecks = TRUE,
                plot = FALSE
                ){
  

    # fix secondary parameters, get implicit parameters
    implicitParameters <- calculateImplicitParameters(
          # gradient scale parameters
          eventChangeScale = eventChangeScale,
          bgGradientValue = bgGradientValue,  
          fullGradientRange = fullGradientRange,
    
          # event to sample scale ratio
          eventSampleWidthRatio = eventSampleWidthRatio,
          # initial secondary parameters
          sampleWidth = sampleWidth,
          eventDuration = eventDuration,
          sedRatePerTimestep = sedRatePerTimestep,
        
          maxSampleTimeStep = maxSampleTimeStep,
          minSampleTimeStep = minSampleTimeStep,
    
          # additional primary paramters
          samplingCompleteness = samplingCompleteness,
          transitionDurationRatio = transitionDurationRatio,
          bioturbDepthRatio = bioturbDepthRatio
          )
    
    # replace secondary parameters
    #eventSampleWidthRatio <- implicitParameters$eventSampleWidthRatio
    #sampleWidth <- implicitParameters$sampleWidth
    #eventDuration <- implicitParameters$eventDuration
    #sedRatePerTimestep <- implicitParameters$sedRatePerTimestep
        
    # Set Up the Pattern of Simulated Gradient Change Over Time
    simGradientChangeOut <- setupSimulatedGradientChange(
          nEvents = nEvents,
          fullGradientRange = fullGradientRange,
          peakGradientValue = implicitParameters$peakGradientValue, 
          bgGradientValue = bgGradientValue,
          bgDurationRange = implicitParameters$bgDurationRange,
          transitionDuration = implicitParameters$transitionDuration,
          eventDuration = implicitParameters$eventDuration,
          halfGradientOnly = halfGradientOnly,
          includeInitialSetUpPhases = includeInitialBackgroundPhase,
          plot = plot
          )
    
    #approxGradientSeriesFunction <- approxGradientSeriesFunction
    
    maxTime <- max(simGradientChangeOut$simGradient$time)
    
    # pulling variables I don't actually need...
    #eventPhaseStartTimes <- simGradientChangeOut$eventPhaseStartTimes
    #backgroundStartEnd <- simGradientChangeOut$backgroundStartEnd
    # What are the start times for each seperate event phase?
      # eventPhaseStartTimes
    # And start-end for the background interval?
      # backgroundStartEnd
    
    ###########################################################
    # Figure out which samples will be taken and at what core depths / outcrop heights
    # First, make a matrix with age, depth, gradient values
    simTimeVar <- data.frame(
        timestep = 1:maxTime,
        # reverse core depth so oldest time is at bottom (biggest depths)
          # 07-08-21: subtract one to start from 0
        sedColumnDepth = ((1:maxTime) - 1) * implicitParameters$sedRatePerTimestep,
        gradientValue = simGradientChangeOut$
            approxGradientSeriesFunction(1:maxTime)
        )
    
    # Simulate actual specimen abundances for each time step
    # using KDEs and occurrence probabilities from empirical data
    timestepAbundances <- getTimestepAbundances(
        kdeRescaled = kdeRescaled,
        specimensPerTimestep = specimensPerTimestep,
        probSpeciesOccur = probSpeciesOccur,
        gradientValues = simTimeVar$gradientValue
        )
    colnames(timestepAbundances) <- names(kdeRescaled)
    
    # Running the Simulation
    fossilSeries <- sampleFossilAssemblageSeries(    
        bioturbIntensity = bioturbIntensity, 
        bioturbZoneDepth = implicitParameters$bioturbZoneDepth,  
        distBetweenSamples = implicitParameters$distBetweenSamples, 
        sampleWidth = implicitParameters$sampleWidth, 
        simTimeVar = simTimeVar,
        timestepAbundances = timestepAbundances,
        nSpecimens = nSpecimens
        )

    # checks for debugging    
    # total number sampled for each species
    #colSums(abundanceTable)[1:10]
    # total number of specimens in each sediment sample
    #rowSums(abundanceTable)[1:10]

    # Check abundances and make sure there is: 
      # (a) sufficient sample in each picked sample, 
      # (b) not an excessive sample in each picked sample, 
      # (c) samples aren't over-dominated by a single taxon. 
    # The first two issues shouldn't arise as problems 
      # when we use a fixed, non-stochastic sample size.
    
    #if(runChecks){
    #    fossilSeries <- checkFossilSeriesAbundanceTable(
    #        fossilSeries = fossilSeries, 
    #        minPickedSampleSize = minPickedSampleSize, 
    #        maxPickedSampleSize = maxPickedSampleSize,
    #        maxDominance = maxDominance
    #        )
    #    }

    ###################################################################
    # Prepping Simulation Data for Post- Simulation Analysis
    sampleProperties <- getSampleProperties(
        simTimeVar = fossilSeries$simTimeVar, 
        fossilSeries = fossilSeries,
        eventStartEndTimes = simGradientChangeOut$eventStartEndTimes,
        initialBackgroundIntervalIncluded = includeInitialBackgroundPhase,
        backgroundStartEnd = simGradientChangeOut$backgroundStartEnd
        )
    
    # check event sampling
    if(samplingCompleteness == 1){
        # test if any unsampled events
            uniqueEvents <- unique(sampleProperties$eventID[!is.na(sampleProperties$eventID)])
            if(length(uniqueEvents) != nEvents){
                stop("Fully sampled record but not all events sampled -- that's impossible Dave!")
               }
        }

    
    # get sample DCA1 scores
    ecologyOutList <- quantifyCommunityEcology(
        origAbundData = origAbundData,
        fossilSeries = fossilSeries,
        useTransformedRelAbundance = useTransformedRelAbundance,
        projectIntoOrigDCA = projectIntoOrigDCA,
        powerRootTransform = powerRootTransform, 
        singularDCA = singularDCA,
        inclusiveDCA = inclusiveDCA,
        rawDCA = rawDCA
        )
    
    # add DCA1 scores to sample properties
    if(singularDCA){
        sampleProperties$scoreDCA1_singular <- ecologyOutList$scoreDCA1_singular 
        }
    if(inclusiveDCA){
        sampleProperties$scoreDCA1_inclusive <- ecologyOutList$scoreDCA1_inclusive 
        }
    if(rawDCA){
        sampleProperties$scoreDCA1_raw <- ecologyOutList$scoreDCA1_raw 
        }
        
    ##########################################################################
    outList <- list(
        implicitParameters = implicitParameters,
        simGradientChangeOut = simGradientChangeOut,
        maxTime = maxTime,
        simTimeVar = simTimeVar,
        fossilSeries = fossilSeries,
        ecology = ecologyOutList,
        sampleProperties = sampleProperties
        )
    
    if(plot){
        plotFossilAssemblageSeriesDCA(
            simTimeVar = outList$simTimeVar, 
            scoreDCA1 = outList$ecology$scoreDCA1_singular,
            sampleAge = outList$sampleProperties$sampleMidAge
            )
        }
    
    return(outList)
    }
