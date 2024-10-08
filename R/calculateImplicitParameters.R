#' Calculate Implicit Parameters for Modeling Time-Series of Fossil Assemblages
#' 
#' Given a sufficient set of parameters for simulating fossil assemblages 
#' as a time-series, this function calculates the full set of parameters
#' necessary for running each component model.

#  using a particular while some parameterizations

#' @details 
#' Under the models considered in \code{paleoAM}, some parameterizations 
#' may be equivalent, even though the a particular analysis might 
#' be better simulated using a particular set of parameters. 
#' Allowing various different parameterizations is a useful 
#' generalization, but requires translating those equivalent parameters 
#' from one set to another (e.g. specifying parameters A, C & D, 
#' but running a simulation that requires parameters A, B & C).
#' 
#' This function mainly exists to calculate unspecified parameters 
#' for \code{\link{simulateFossilAssemblageSeries}} and also to 
#' identify conflicting parameter specifications.

#' @inheritParams setupSimulatedGradientChange

#' @param eventChangeScale A value indicating the amount relative to 
#' the background value (\code{bgGradientValue}) and the maximum 
#' possible change as indicated by \code{fullGradientRange} 
#' (in other words, simulated change must be within observed gradient, 
#' so \code{eventChangeScale} is a proportional multiplier of 
#' the total possible change).

#' @param fullGradientRange A vector of two values giving the minimum 
#' and maximum gradient values observed in the empirical data.

#' @param eventSampleWidthRatio How long should an event be relative 
#' to the amount of time (or sediment) captured within a sedimentary sample? 
#' This parameter is used for simulating event duration, 
#' sample width and sedimentation rate where any two of these three 
#' are defined and the third is not defined. 
#' This value is referred to as \emph{Resolution Potential} 
#' in Belanger & Bapst (2023).

#' @param sampleWidth The 'width' of a sample relative to core depth 
#' or outcrop height, usually given in linear units (usually centimeters). 
#' For taking sediment samples from a core, this is straightforward 
#' (how thick is each sediment sample taken?) but for outcrops this 
#' may be more difficult to determine 
#' (what is the thickness of a horizon in a shale unit?).

#' @param sedRatePerTimestep The rate of sedimentation, given as a 
#' ratio of sediment thickness (given in linear dimensions, 
#' in the same units as \code{sampleWidth}), over time 
#' (given in the same time units as \code{eventDuration}. 

#' @param maxSampleTimeStep The maximum number of individual time-steps 
#' used for simulating a sample.

#' @param minSampleTimeStep The minimum number of individual time-steps 
#' used for simulating a sample.

#' @param samplingCompleteness The relative completeness of stratigraphic 
#' sampling. For example, if two-centimeter wide samples of sediment are 
#' taken from a sediment core, every ten centimeters, then the 
#' \code{samplingCompleteness} is two over 10, or that
#' \code{samplingCompleteness = 1/5}. A simulation with a sampling 
#' completeness of 1 would be comparable to exhaustively sampling a 
#' core that recorded no gaps in sedimentation over its history. 
#' Rocky outcrops are more complicated, as fossil-bearing horizons 
#' may be relatively thin compared to the thickness of the section, 
#' such that outcrop-based fossil records should be simulated as 
#' having \emph{very low} \code{samplingCompleteness}.
#' 

#' @param transitionDurationRatio The ratio of how long the transition 
#' between peak and background intervals should be, relative to the 
#' length of the peak 'event' duration (\code{eventDuration}). 
#' The longer this transition interval, the more chances of 
#' an assemblage being sampled that represents transitional gradient values.

#' @param bioturbDepthRatio The ratio of the sediment depth to which 
#' bioturbation occurs, made relative to the width of a 
#' sediment sample (\code{sampleWidth}). 
#' A \code{sampleWidth} of 3 cm and a \code{biotubDepthRatio} 
#' of 5 implies a bioturbation depth of 15 cm (\code{3 * 5}).  
#' Bioturbation depth varies considerably in the modern ocean, but is 
#' often the depth of active bioturbation is about 10 cm, such that the 
#' top ten centimeters of sediment 
#' (and the organic remains in those ten centimeters of sediment) 
#' are being regularly moved up and down by organism activity. 
#' For the purposes of this model, a bioturbation zone depth of 
#' 10 centimeters means that sampling a centimeter of sediment 
#' at location X, the apparent fossil assemblage that would be 
#' recovered is just as likely to include specimens that were 
#' deposited five centimeters away as those deposited at location X.

#' @return
#' Returns a list giving the full set of parameters necessary for running \code{\link{simulateFossilAssemblageSeries}}.

# @aliases

#' @seealso
#' \code{\link{simulateFossilAssemblageSeries}}

#' @references
#' Belanger, Christina L., and David W. Bapst. 2023.
#' "Simulating our ability to accurately detect abrupt 
#' changes in assemblage-based paleoenvironmental proxies." 
#' Palaeontologia Electronica 26 (2), 1-32

# @examples


#' @name calculateImplicitParameters
#' @rdname calculateImplicitParameters
#' @export
calculateImplicitParameters <- function(
        # gradient scale parameters
        eventChangeScale,
        bgGradientValue,  
        fullGradientRange,
  
        # event to sample scale ratio
        eventSampleWidthRatio = NULL,
        # initial secondary parameters
        sampleWidth = NULL,
        eventDuration = NULL,
        sedRatePerTimestep = NULL,
        # number of time-step assemblages
            # to simulate per sample
        maxSampleTimeStep = 500,
        minSampleTimeStep = 3,
  
        # additional primary paramters
        samplingCompleteness,
        transitionDurationRatio,    
        bioturbDepthRatio
        
        ){

    # check sampling completeness
    if(!(samplingCompleteness > 0)){
        stop("samplingCompleteness must be more than zero")
        }
    
    
    # fixing secondary parameters
    
    # To calculate the necessary parameters, we need to set *THREE* of *FOUR* arbitrary 'nuisance' parameters: 
      # **sample width** (in cm),
      # **duration** that an event 'plateaus' at some peak value (in time-steps) 
      # **sedimentation rate** (in cm per time-steps)
      # **sample/peak width ratio** (the ratio of the length of the event relative to the sample width
    # Combining the three known variables composes a closed system, 
      # and the fourth unknown can be calculated from the other three variables.
    # Ultimately the choice of values for the above do not matter to many analyses, 
      # as phenomena should scale with our key parameter, 
      # the relationship between sample width and peak width.
    
    nAbsentSecParam <- sum(c(
        is.null(eventSampleWidthRatio), 
        is.null(sampleWidth), 
        is.null(eventDuration), 
        is.null(sedRatePerTimestep)
        ))
    if(nAbsentSecParam < 1){
        stop(paste0("eventSampleWidthRatio and the three", 
            " related secondary parameters (sampleWidth, eventDuration and sedRatePerTimestep)\n",
            " cannot *all* be given input values, as any three constrain the fourth."))
        }
    if(nAbsentSecParam > 1){
        stop(paste0("Three, and only three, of the variables eventSampleWidthRatio\n", 
            " and related secondary parameters (sampleWidth, eventDuration and sedRatePerTimestep)\n",
            " must be given input values, so to constrain the fourth."))
        }

    if(is.null(eventSampleWidthRatio)){
        eventSampleWidthRatio <- sedRatePerTimestep * eventDuration / sampleWidth
        }
        
    if(is.null(sampleWidth)){
        sampleWidth <- sedRatePerTimestep * eventDuration / eventSampleWidthRatio
        }
    
    if(is.null(eventDuration)){
        eventDuration <- sampleWidth * eventSampleWidthRatio / sedRatePerTimestep 
        }
    
    if(is.null(sedRatePerTimestep)){
        # Given the above, we can calculate the effective sedimentation rate as:
        sedRatePerTimestep <- sampleWidth * eventSampleWidthRatio / eventDuration 
        }
    
    # adjust all four of above so we don't simulate too many time-steps per sample
    #multPar <- maxSampleTimeStep / (sampleWidth/sedRatePerTimestep) 
    # adjust sampling rate or sampleWidth? Hmmmmm. Both?
    #sampleWidth <- sampleWidth * multPar
    #sedRatePerTimestep <- multPar/sedRatePerTimestep
    #eventSampleWidthRatio <- 
    expStepsPerSample <- (sampleWidth/sedRatePerTimestep)
    if(maxSampleTimeStep < expStepsPerSample){
        stop("More time steps expected in a sample than maxSampleTimeStep -- may be too computationally intensive")
        }
    if(minSampleTimeStep > expStepsPerSample){
        stop("Fewer time steps expected in a sample than minSampleTimeStep -- may be too computationally intensive")
        }
    
    nAbsentSecParam <- sum(c(
        is.null(eventSampleWidthRatio), 
        is.null(sampleWidth), 
        is.null(eventDuration), 
        is.null(sedRatePerTimestep)
        ))
    
    if(nAbsentSecParam > 0){
        stop("Unexpected NULL parameter for sed rate still NULL ! Please investigate.")
        }

    # Implicit Parameters (Not Directly Set)
    
    ## Calculate Actual Durations of Samples, Events, Etc
    
    # The implicit duration of samples is given by:
    
    sampleDuration <- eventDuration / eventSampleWidthRatio
    
    # As described above, we use *sampling completeness*, 
      # the ratio of the sample width to the total sampling interval. 
    # Thus, we can calculate the distance between samples as:
    
    # calculate the total length of the sampling interval width
    samplingIntervalWidth <- sampleWidth/samplingCompleteness
    
    # calculate the distance between samples as the difference between
      # sampling interval width and the sample width
    distBetweenSamples <- samplingIntervalWidth - sampleWidth
    
    # The depth of the bioturbation zone is given by 
      # the product of the bioturbation depth ratio and the sample width.
    
    # Size of region that sediment particle mixture occurs in (in centimeters)
    bioturbZoneDepth <- bioturbDepthRatio * sampleWidth
    
    # how long should it take for values to transition from background to peak values?
    transitionDuration <- transitionDurationRatio * eventDuration 
    
    # if small, then transition is very fast, to avoid having much sampling of communities
      # in between background and peak
    
    # The duration of events in simulation time-units 
      # (time-steps, sometimes referred to as 'years') 
      # is defined by the secondary parameter 
      # for the duration of each event's plateau / peak (`eventDuration`). 
    # From this and the transition duration ratio, 
      # we can calculate the duration of each transition interval, 
      # leading into and out of each event.
        
    # How long should the record be? Well, we need 
      # somewhat stable background phases at the beginning, the end, and between pulsed spikes. 
    # However, if we make these have the same identical duration, 
      # then we get artificial patterns of which peaks we sample 
      # that are identical to our patterns of how we took sample from the core / outcrop. 
    # That makes sense but also is completely weird-looking and hard to explain. 
    
    # So we'll add a small stochastic element to these simulations by pulling 
      # those waiting times from a uniform distribution. 
    # Thus, background intervals are of variable duration.
    
    # set duration of background times as a uniform distribution
      # to avoid artificial patterns of hits and misses for spikes
    # Constrain so buffer around events is always at as large 
      # as the effective duration represented by a sample's width
    baseDurationBG <- max(
        c(eventDuration, sampleDuration, bioturbZoneDepth, distBetweenSamples)
        ) * 1.1

    minBgDuration <- baseDurationBG * 2
    maxBgDuration <- baseDurationBG * 3
    bgDurationRange <- c(minBgDuration, maxBgDuration)
    
    # At `r sedRatePerTimestep` cm/yr sedimentation rate, 
      # that would mean that there is (on average) 
      # `r (minBgDuration + maxBgDuration)/2*sedRatePerTimestep` cm in a background interval.
    
    # The sampling resolution of each event is, much like the sedimentation rate, 
      # an emergent property of the model when other parameters are defined,
      # depending on both the effect size of sample widths, 
      # and the amount of unsampled sediment between samples. 
    # Given our parameterization, it makes sense to talk about it in terms of 
      # how many sampling instances are expected to fall within an event duration. 
    # This is a unit-less value, simply the product of the ratio of sampling width to event width
      # and the ratio of sample width to sampling interval width.
    
    samplingEventResolution <- samplingCompleteness / eventSampleWidthRatio
    #samplingEventResolution

    ## background and event peak values
    
    fullGradientDiff <- abs(diff(fullGradientRange))
    
    # peak gradient values at times of simulated spikes
    peakGradientValue <- (eventChangeScale * fullGradientDiff) + bgGradientValue
    
    if(bgGradientValue < fullGradientRange[1]){
      stop("background gradient value is less than min gradient value")
      }
    if(bgGradientValue > fullGradientRange[2]){
      stop("background gradient value is more than max gradient value")
      }
    if(peakGradientValue < fullGradientRange[1]){
      stop("peak gradient value is less than min gradient value")
      }
    if(peakGradientValue > fullGradientRange[2]){
      stop("peak gradient value is more than max gradient value")
      }

    ##############
    outputList <- list(
        # output fixed ESWR and secondary parameters
        eventSampleWidthRatio = eventSampleWidthRatio,
        sampleWidth = sampleWidth,
        eventDuration = eventDuration,
        sedRatePerTimestep = sedRatePerTimestep,
        # output implicit parameters
        peakGradientValue = peakGradientValue,
        sampleDuration = sampleDuration,
        transitionDuration = transitionDuration,
        # samplingIntervalWidth = samplingIntervalWidth,
        bioturbZoneDepth = bioturbZoneDepth,
        bgDurationRange = bgDurationRange,
        distBetweenSamples = distBetweenSamples,
        samplingEventResolution = samplingEventResolution
        )
    
    return(outputList)
    }
