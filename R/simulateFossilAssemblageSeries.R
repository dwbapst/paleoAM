#' Simulate Time-Series of Successive Fossil Assemblages
#'
#' Given a set of parameters and models describing species abundance,
#' stochastically models changes in an underlying biotic gradient and simulates
#' ecological change and a sequence of samples representing change in recovered
#' fossil assemblages over that interval, 
#' including estimating the recovered gradient. 
#  via \code{quantifyCommunityEcology}.

#' @details
#' Different parameterizations may be given as input, allowing different parameters to be unspecified.
#' Missing paramters are then calculated from the specified ones using \code{\link{calculateImplicitParameters}}.

#' @inheritParams setupSimulatedGradientChange 
#' @inheritParams calculateImplicitParameters 
#' @inheritParams sampleFossilSeries
#' @inheritParams getTimestepAbundances
#' @inheritParams getSampleDCA 
# @inheritParams quantifyCommunityEcology

#' @param specimensPerTimestep The number of specimens returned in a given time-step by \code{getTimestepAbundances}, usually set to an unrealistically high number to represent the true 'unsampled' fossil assemblage. Default is 10000.

#' @param plot Should the simulated time-series of fossil assemblages be shown as a sequence of generating and recovered gradient values against time? Default is \code{FALSE}.
 
#' @param thinOutput Should the output be thinned to just the sample properties and intrinsic variables? Default is FALSE.

#' @return
#' Returns a list, which by default has seven components: 
#' \code{implicitParameters}, the full list of parameters used for generating the simulated data; 
#' \code{simGradientChangeOut}, the simulated time-series of gradient change output by \code{setupSimulatedGradientChange};
#' \code{maxTime}, the total duration of the entire simulated time-series from start to end;
#' \code{simTimeVar}, a data frame specifying time-steps, sedimentary depth and environmental gradient values for simulating a time-series of sampled fossil assemblages, used as input in \code{\link{sampleFossilSeries}};
#' \code{fossilSeries}, a list containing the simulated time-series of sampled fossil assemblages from \code{\link{sampleFossilSeries}},
#' \code{ecology}, the recovered ecological variables for each simulated sample, 
#'      as returned by internal function \code{quantifyCommunityEcology},
#' and \code{sampleProperties}, a list containing a number of variables specific to individual .
#' 
#' If \code{thinList = TRUE} is used, then the output list
#'  contains only two components: 
#'  \code{sampleProperties} and \code{implicitParameters}.
#'  The \code{implicitParameters} component is the same as in the full output,
#'  but the \code{sampleProperties} component only contains information on when 
#'  (in both time and sedimentary depth) a given sample is located in the 
#'  simulated time-series, and the variable \code{scoreDCA1_recovered}. 

# @aliases

#' @seealso
#' \code{\link{calculateImplicitParameters}}

#' @references
#' Belanger, Christina L., and David W. Bapst. 
#' "Simulating our ability to accurately detect abrupt 
#' changes in assemblage-based paleoenvironmental proxies." (2023): 1-32.
#' https://doi.org/10.1073/pnas.1602102113

#' @examples
#' # an example with Gulf of Alaska data
#' \donttest{
#' # load data
#' data(gulfOfAlaska)
#' 
#' alaskaKDEs <- getSpeciesSpecificRescaledKDE(
#'     gradientOrigDCA = DCA1_GOA, 
#'     origAbundData = abundData_GOA, 
#'     abundanceFloorRatio = 0.5, 
#'     nBreaksGradientHist = 20, 
#'     modeledSiteAbundance = 10000
#'     )
#'     
#' alaskaProbOccur <- getProbOccViaPresAbs(
#'    gradientOrigDCA = DCA1_GOA, 
#'    origAbundData = abundData_GOA
#'    )
#' 
#' # Run the simulation of fossil assemblages
#'     # simulateFossilAssemblageSeries has lots of arguments...
#'     # below they are broken up into groups, seperate by #
#'     # matches scenarios from fig 13 of Belanger & Bapst
#'     
#' fossilSeriesOut <- simulateFossilAssemblageSeries(
#'       # inputs
#'       kdeRescaled = alaskaKDEs,
#'       probSpeciesOccur = alaskaProbOccur,
#'       origAbundData = abundData_GOA,
#'       fullGradientRange = c(min(DCA1_GOA), max(DCA1_GOA)),
#'       
#'       # let's make it relatively mild event 
#'         # with a long transition 
#'       eventChangeScale = 0.5,
#'       bgGradientValue = -1,
#'       transitionDurationRatio = 1,
#'        
#'       # don't need to define eventSampleWidthRatio 
#'         # - only need to define three of eventSampleWidthRatio, 
#'         # sampleWidth, eventDuration, sedRatePerTimestep
#'       sampleWidth = 3,
#'       eventDuration = 100, 
#'       sedRatePerTimestep = 0.1,
#'       
#'       # sample every third sample-width worth of core
#'       samplingCompleteness = 1/3,
#'       # no bioturbation 
#'       bioturbDepthRatio = 0,
#'       bioturbIntensity = 0,
#'            
#'       nEvents = 1,
#'       nSpecimens = 100,
#'       # let's plot it
#'       plot = TRUE
#'       )
#' }
#' 

#' @name simulateFossilAssemblageSeries
#' @rdname simulateFossilAssemblageSeries
#' @export
simulateFossilAssemblageSeries <- function(
                kdeRescaled,
                probSpeciesOccur,
                origAbundData,
                eventChangeScale,
                bgGradientValue,
                fullGradientRange,
                eventSampleWidthRatio = NULL,
                sampleWidth = NULL,
                eventDuration = NULL,
                sedRatePerTimestep = NULL,
                samplingCompleteness,
                transitionDurationRatio,
                bioturbDepthRatio,
                bioturbIntensity,     
                nEvents,
                nSpecimens,
                #
                specimensPerTimestep = 10000,
                halfGradientOnly = FALSE,
                useTransformedRelAbundance = TRUE,
                projectIntoOrigDCA = TRUE,
                powerRootTransform = 1,  
                maxSampleTimeStep = 500,
                minSampleTimeStep = 3,
                includeInitialBackgroundPhase = FALSE,
                # runChecks = TRUE,
                plot = FALSE,
                thinOutput = FALSE
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
          peakGradientValue = implicitParameters$peakGradientValue, 
          bgGradientValue = bgGradientValue,
          bgDurationRange = implicitParameters$bgDurationRange,
          transitionDuration = implicitParameters$transitionDuration,
          eventDuration = implicitParameters$eventDuration,
          halfGradientOnly = halfGradientOnly,
          includeInitialBackgroundPhase = includeInitialBackgroundPhase,
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
    fossilSeries <- sampleFossilSeries(    
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
    #        maxDominance = maxDominance,
    #        speciesNames = speciesNames
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
        powerRootTransform = powerRootTransform
        )
    
    # add DCA1 scores to sample properties
        # we now only do singular / projection method for DCA
    
    sampleProperties$scoreDCA1_recovered <- ecologyOutList$scoreDCA1_recovered 
        
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
            gradientRecovered = outList$ecology$scoreDCA1_recovered,
            sampleAge = outList$sampleProperties$sampleMidAge
            )
        }
    
    thinList <- function(input){
        out <- list(sampleProperties = list(), implicitParameters = list())
        out$sampleProperties$sampleSedColumnDepth_start <- input$sampleProperties$sampleSedColumnDepth_start
        out$sampleProperties$sampleSedColumnDepth_end <- input$sampleProperties$sampleSedColumnDepth_end
        out$sampleProperties$sampleInterval_start  <- input$sampleProperties$sampleInterval_start 
        out$sampleProperties$sampleInterval_end <- input$sampleProperties$sampleInterval_end   
        out$sampleProperties$sampleMidAge <- input$sampleProperties$sampleMidAge
        out$sampleProperties$scoreDCA1_recovered <- input$sampleProperties$scoreDCA1_recovered
        out$implicitParameters <- input$implicitParameters
        return(out)
        }

    if(thinOutput){
        outList <- thinList(outList)
        }
        
    return(outList)
    }
