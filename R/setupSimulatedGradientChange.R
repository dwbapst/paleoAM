#' Create a Stochastic Time-Series of Gradient Change For Use in Simulating Assemblage Change
#'
#' Given a series of inputs, simulates a sequence of gradient change 
#' against time for use in testing how environmental change alters 
#' the recovered sequence of fossil assemblages.

#' @details 
#' This function is rather complicated and was written at a time 
#' when it was envisioned that simulations would involve time series 
#' of many repeated events with varying background intervals 
#' between them, rather than simulated sequences having only one event. 
#' In practice, use of paleoAM has tended to find the latter to be more useful.

#' @param nEvents Number of events to occur in a simulated sequence 
#' of gradient change.

#' @param peakGradientValue The gradient value at the 'peak' for 
#' an event that represents an excursion on that environmental gradient.

#' @param bgGradientValue The gradient value expected during 
#' background intervals during which no notable excursion is 
#' occurring on that environmental gradient.

#' @param bgDurationRange A vector of two values, representing the 
#' minimum and maximum duration (in time units) of a 
#' background interval between successive events.

#' @param transitionDuration How long the transition between peak and 
#' background intervals should be. The longer this interval, 
#' the more chances of an assemblage being sampled that 
#' represents transitional gradient values.

#' @param eventDuration The duration (in time-units) of a 
#' simulated event during which the environmental gradient 
#' is at an excursion 'peak' level.

#' @param halfGradientOnly Whether to simulate only half of 
#' a background-event sequence, either beginning or terminating 
#' the simulation at the peak value. 
#' Only a single event can be simulated, so \code{nEvents} must be 1. 
#' The default is \code{FALSE} which signals to not simulated a half-gradient.

#' @param includeInitialBackgroundPhase A logical indicating whether 
#' to include a lengthy background phase, for use in calibrating a simulation. 
#' This function is mainly for diagnostic purposes 
#' and may be removed in future updates.

#' @param plot Should the simulated gradient be shown as a plot?

#' @return 
#' A list with five components: 
#' \code{simGradient}, a data frame giving the change in gradient values over time; 
#' \code{approxGradientSeriesFunction}, the simulated gradient given as an interpolated function;
#' \code{eventStartEndTimes}, a vector of when each event and its preceding transition begin in time-units;
#' \code{eventPhaseStartTimes}, a vector of when each new event phase (at the peak gradient value) begin in time-units;
#' and \code{backgroundStartEnd}, a value indicating the time-step when the beginning background interval ends.
#' 

# @examples


#' @name setupSimulatedGradientChange
#' @rdname setupSimulatedGradientChange
#' @export
setupSimulatedGradientChange <- function(
                nEvents,
                peakGradientValue, 
                bgGradientValue,
                bgDurationRange,
                transitionDuration,
                eventDuration,
                halfGradientOnly = FALSE,
                includeInitialBackgroundPhase = TRUE,
                plot = FALSE
                ){
    
    # Set Up the Pattern of Simulated Gradient Change Over Time
    
    # This will effectively be a time-series of gradient values, 
    # with irregular intervals between dates.
    
    if(halfGradientOnly == FALSE){
        if(nEvents != 1){
            stop("halfGradientOnly options can only be used if simulating a single event")
            }
        }
    
    # Background intervals are of uneven duration and thus will be obtained using runif
    bgDuration <- function(n=1, min = bgDurationRange[1], max = bgDurationRange[2]){
        stats::runif(n, min = min, max = max)
        }
    
    if(includeInitialBackgroundPhase){
    
        # The initial segment of the gradient timeline is used to calibrate later results, 
        # with both an instantaneous and a slow transition between two minimum and maximum values. 
        # The building of the timeline is iterative and needs to be done step-by-step 
        # as the length of some intervals depend on the length of previous interval durations.
        
        # first simulate a plateau at the peak, and at background
        simGradientTime <- c(0, bgDuration()*3)
        # fast transition from peak to background
        simGradientTime <- c(simGradientTime, max(simGradientTime) * 1.001)
        simGradientTime <- c(simGradientTime, max(simGradientTime) + bgDuration()*3)
        # simulate a gradual slope-up and slope down
        simGradientTime <- c(simGradientTime, max(simGradientTime) + bgDuration()*3)
        simGradientTime <- c(simGradientTime, max(simGradientTime) + bgDuration()*3)
        simGradientTime <- c(simGradientTime, max(simGradientTime) + bgDuration()*3)
        
        # long background interval
            # record start and end for calculating a background envelope 
            # but buffer before and after recording
        # start buffer
        simGradientTime <- c(simGradientTime, max(simGradientTime) + bgDuration())
        # record time
        backgroundStartEnd <- max(simGradientTime)
        simGradientTime <- c(simGradientTime, max(simGradientTime) + bgDuration()*20)
        # record time
        backgroundStartEnd <- c(backgroundStartEnd, max(simGradientTime))
        # end buffer
        simGradientTime <- c(simGradientTime, max(simGradientTime) + bgDuration())
        
        # Assign gradient values to the gradient time-series.
        
        # get gradient values at inflection points
        # peak, peak, background, background, peak, peak, background, background
        simGradientValue <- c(
            peakGradientValue, peakGradientValue,
            bgGradientValue, bgGradientValue,
            peakGradientValue, peakGradientValue,
            bgGradientValue, bgGradientValue,
            bgGradientValue, bgGradientValue
            )
    }else{
        # just start simulation with a single background interval
        simGradientTime <- c(0, bgDuration() * 2)
        simGradientValue <- c(bgGradientValue, bgGradientValue)
        backgroundStartEnd <- c(NA, NA)
        }
    
    # record when spike simulation begins
    eventPhaseStartTimes <- numeric(length = nEvents)
    eventStartEndTimes <- matrix(NA, nEvents, 2)
    
    # Now add in the spikes at the specified peak and background height.
    for(i in 1:nEvents){
        lastTime <- max(simGradientTime)
        eventPhaseStartTimes[i] <- lastTime
        newGradientTime <- lastTime + bgDuration()
        lastTime <- newGradientTime
        
        if(halfGradientOnly == FALSE){
            newGradientTime <- c(newGradientTime,
                lastTime + transitionDuration,
                lastTime + transitionDuration + eventDuration,
                lastTime + transitionDuration + eventDuration + transitionDuration,
                lastTime + transitionDuration + eventDuration + transitionDuration + bgDuration()
                )
            # check
            if(any(is.na(newGradientTime[2:3]))){
                stop("NAs in newGradientTime in setupSimulatedGradientChange")
                }
            # 
            eventStartEndTimes[i,] <- newGradientTime[2:3]
            simGradientTime <- c(simGradientTime, newGradientTime)
            #
            newGradientValue <- c(bgGradientValue,
                peakGradientValue, peakGradientValue,
                bgGradientValue, bgGradientValue)
            simGradientValue <- c(simGradientValue, newGradientValue)
            
            }else{
                if(halfGradientOnly == "riseOnly"){
                    
                    newGradientTime <- rev(c(
                        newGradientTime,
                        lastTime + transitionDuration,
                        lastTime + transitionDuration + eventDuration
                        ))
                    
                    # check
                    if(any(is.na(newGradientTime[2:3]))){
                        stop("NAs in newGradientTime in setupSimulatedGradientChange")
                        }
                    # 
                    eventStartEndTimes[i,] <- newGradientTime[2:3]
                    simGradientTime <- c(newGradientTime, simGradientTime)
                    #
                    newGradientValue <- rev(c(bgGradientValue,
                        peakGradientValue, peakGradientValue))
                    
                    simGradientValue <- c(newGradientValue, simGradientValue)
                    }
                }
        }
    
    # get total gradient time -> maxTime
    #maxTime <- max(simGradientTime)
    #meanBgDuration <- (minBgDuration + maxBgDuration)/2
    # Now we will use `approx` to build a function that will give us the gradient value
    # for each point in time in the simulation.
    
    if(length(simGradientTime) != length(simGradientValue)){
        stop(paste0("simGradientTime length (", length(simGradientTime), 
             ") is not equal to simGradientValue length (",
             length(simGradientValue),")"
             ))
        }

    # checks
    if(any(is.na(eventStartEndTimes))){
        stop("NAs in eventStartEndTimes created by setupSimulatedGradientChange")
        }
    if(any(is.na(eventPhaseStartTimes))){
        stop("NAs in eventPhaseStartTimes created by setupSimulatedGradientChange")
        }    
    #if(any(is.na(backgroundStartEnd))){
    #    stop("NAs in backgroundStartEnd created by setupSimulatedGradientChange")
    #    }
    if(any(is.na(simGradientTime))){
        stop("NAs in simGradientTime created by setupSimulatedGradientChange")
        }       

    # 07-26-21
    # invert time scale so timesteps increases with depth
        # like in a real sedimentary core, you nut job
    eventStartEndTimes   <- max(simGradientTime) - eventStartEndTimes
    eventPhaseStartTimes <- max(simGradientTime) - eventPhaseStartTimes
    backgroundStartEnd   <- max(simGradientTime) - backgroundStartEnd
    simGradientTime      <- max(simGradientTime) - simGradientTime
    
    # build approximation function for gradient values
    approxGradientSeriesFunction <- function(time){
        stats::approx(x = simGradientTime, 
               y = simGradientValue, 
               xout = time)$y
        }
    
    if(plot){
        plot(1:max(simGradientTime), 
             approxGradientSeriesFunction(1:max(simGradientTime)), 
             xlim = c(max(simGradientTime), 0),
             type = "l",
             main = "Simulated Change in Gradient Values over Time",
             xlab = "Simulation Time", 
             ylab = "Simulated Gradient Values"
            )
        }
    
    output <- list(
        simGradient = data.frame(time = simGradientTime, value = simGradientValue),
        approxGradientSeriesFunction = approxGradientSeriesFunction,
        eventStartEndTimes = eventStartEndTimes,
        eventPhaseStartTimes = eventPhaseStartTimes,
        backgroundStartEnd = backgroundStartEnd
        )
    
    return(output)
    }