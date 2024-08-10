#' @details

#' @inheritParams

#' @param nEvents 

#' @param fullGradientRange 

#' @param peakGradientValue 

#' @param bgGradientValue 

#' @param bgDurationRange 

#' @param transitionDuration 

#' @param eventDuration 

#' @param halfGradientOnly 

#' @param includeInitialBackgroundPhase 

#' @param plot 

#' @return

#' @aliases

#' @seealso

#' @references

#' @examples


#'
#' @name
#' @rdname
#' @export


setupSimulatedGradientChange <- function(
                nEvents,
                fullGradientRange,
                peakGradientValue, 
                bgGradientValue,
                bgDurationRange,
                transitionDuration,
                eventDuration,
                halfGradientOnly = "full",
                includeInitialBackgroundPhase = TRUE,
                plot = FALSE
                ){
    
    # Set Up the Pattern of Simulated Gradient Change Over Time
    
    # This will effectively be a time-series of gradient values, 
    # with irregular intervals between dates.
    
    if(halfGradientOnly != "full"){
        if(nEvents != 1){
            stop("halfGradientOnly options can only be used if simulating a single event")
            }
        }
    
    # Background intervals are of uneven duration and thus will be obtained using runif
    bgDuration <- function(n=1, min = bgDurationRange[1], max = bgDurationRange[2]){
        runif(n, min = min, max = max)
        }
    
    if(includeInitialSetUpPhases){
    
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
            fullGradientRange[2], fullGradientRange[2], 
            fullGradientRange[1], fullGradientRange[1],
            fullGradientRange[2], fullGradientRange[2], 
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
        
        if(halfGradientOnly == "full"){
            newGradientTime <- c(newGradientTime,
                lastTime + transitionDuration,
                lastTime + transitionDuration + eventDuration,
                lastTime + transitionDuration + eventDuration + transitionDuration,
                lastTime + transitionDuration + eventDuration + transitionDuration + bgDuration()
                )
            eventStartEndTimes[i,] <- newGradientTime[2:3]
            simGradientTime <- c(simGradientTime, newGradientTime)
            #
            newGradientValue <- c(bgGradientValue,
                peakGradientValue, peakGradientValue,
                bgGradientValue, bgGradientValue)
            simGradientValue <- c(simGradientValue, newGradientValue)
            
            }else{
                if(halfGradientOnly == "riseOnly"){
                    
                    newGradientTime <- rev(c(newGradientTime,
                        lastTime + transitionDuration,
                        lastTime + transitionDuration + eventDuration
                        ))
                
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
    
    # 07-26-21
    # invert time scale so timesteps increases with depth
        # like in a real sedimentary core, you nut job
    eventStartEndTimes <- max(simGradientTime) - eventStartEndTimes
    eventPhaseStartTimes <- max(simGradientTime) - eventPhaseStartTimes
    backgroundStartEnd <- max(simGradientTime) - backgroundStartEnd
    simGradientTime <- max(simGradientTime) - simGradientTime
    
    # build approximation function for gradient values
    approxGradientSeriesFunction <- function(time){
        approx(x = simGradientTime, 
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