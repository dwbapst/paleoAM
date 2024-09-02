#' Calculate Implicit Parameters 

#' @details 

#' @inheritParams setupSimulatedGradientChange

#' @param eventChangeScale A value indicating the amount relative to the background value (\code{bgGradientValue}) and the maximum possible change as indicated by \code{fullGradientRange} (in other words, simulated change must be within observed gradient, so \code{eventChangeScale} is a proportional multiplier of the total possible change).

#' @param fullGradientRange A vector of two values giving the minumum and maximum gradient values observed in the empirical data.

#' @param eventSampleWidthRatio How long should an event be relative to the amount of time (or sediment) captured within a sedimentary sample? This parameter is used for simulating event duration, sample width and sedimentation rate where any two of these three are defined and the third is not defined.

#' @param sampleWidth The 'width' of a sample relative to core depth or outcrop height, usually given in linear units (usually centimeters). For taking sediment samples from a core, this is straightforward (how thick is each sediment sample taken?) but for outcrops this may be more difficult to determine (what is the thickness of a horizon in a shale unit?).

#' @param sedRatePerTimestep The  

#' @param maxSampleTimeStep 

#' @param minSampleTimeStep 

#' @param samplingCompleteness 

#' @param transitionDurationRatio 

#' @param bioturbDepthRatio 

#' @return

#' @aliases

#' @seealso

#' @references

#' @examples


#' @name
#' @rdname
#' @export
calculateImplicitParameters <- function(
        # gradient scale parameters
        eventChangeScale,
        bgGradientValue,  
        fullGradientRange,
  
        # event to sample scale ratio
        eventSampleWidthRatio,
        # initial secondary parameters
        sampleWidth,
        eventDuration,
        sedRatePerTimestep,
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
