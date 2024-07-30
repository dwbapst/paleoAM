#' @details

#' @inheritParams

#' @param

#' @return

#' @aliases

#' @seealso

#' @references

#' @examples

#' @param simTimeVar 
#'
#' @param fossilSeries 
#' @param eventStartEndTimes 
#' @param initialBackgroundIntervalIncluded 
#' @param backgroundStartEnd 
#'
#' @name
#' @rdname
#' @export


getSampleProperties <- function(
            simTimeVar, 
            fossilSeries, 
            eventStartEndTimes,
            initialBackgroundIntervalIncluded,
            backgroundStartEnd
            ){
    
    ###################################################
    # Now, we will get sample age and approximate the true gradient value 
      # for each sample, so we can compare the 'truth' against 
      # what we observe in ordinations applied to the simulated abundance data.
    
    # get ages and the 'true' gradient values
      # approximate the age of each sample based on mid-depth
      # note ages may not be exact timesteps due to rounding issues
    
    # 07-26-21
        # time always runs backwards now, like a real geologic record
        # CB: "Time doesn't run backwards, it runs 'down depth'."
    
    # age mid points
    #sampleMidAge <- approx(
    #    x = simTimeVar$sedColumnDepth, 
    #    y = simTimeVar$timestep, 
    #    xout = fossilSeries$sampleMidDepth
    #    )$y
    
    # age as an interval
    sampleInterval_start <- approx(
        x = simTimeVar$sedColumnDepth, 
        y = simTimeVar$timestep, 
        xout = fossilSeries$sampleIntervals[,1]
        )$y
    
    # and the other interval
    sampleInterval_end <- approx(
        x = simTimeVar$sedColumnDepth, 
        y = simTimeVar$timestep, 
        xout = fossilSeries$sampleIntervals[,2]
        )$y
    
    # combine
    sampleIntervalAges <- cbind(sampleInterval_start, sampleInterval_end)
    
    # check that starts come before ends
    
    if(!all(sampleInterval_start >= sampleInterval_end)){
        stop("starting dates should be larger than end dates")
        }
    
    if(any(is.na(sampleIntervalAges))){
        stop("NAs in sampleIntervalAges")
        }
    
    # also approx the generating gradient values
      # for each sample
    sampleGradientValues <- approx(
        x = simTimeVar$sedColumnDepth, 
        y = simTimeVar$gradientValue, 
        xout = apply(sampleIntervalAges, 1, mean)
        )$y
    
    ##########################################################
    # is a sample in the initial background segment?

    # first... did the simulation even include an initial background segment?
    
    if(initialBackgroundIntervalIncluded){
        isBackgroundSegment <- (sampleIntervalAges[,1] > backgroundStartEnd[2]  &  
            sampleIntervalAges[,2] < backgroundStartEnd[1])
        
        if(length(isBackgroundSegment)<1 | sum(isBackgroundSegment)<1){
            stop("No samples from background interval?")
            }
    }else{
        isBackgroundSegment <- rep(NA, nrow(sampleIntervalAges))
        }
    
    # is a sample overlapping with one of the events? 
    eventID <- apply(sampleIntervalAges, 1, ageMatchFun, 
        eventStartEndTimes = eventStartEndTimes)   

    if(all(is.na(eventID))){
        stop("No samples found during events. Something very bad happened.")
        }
    
    if(any(c(length(sampleInterval_start), 
             length(sampleInterval_end), 
             nrow(fossilSeries$bioturbIntervals), 
             length(sampleGradientValues), 
             length(isBackgroundSegment)
             ) != length(eventID))){
                 stop("sample-wise variables are not identical length")
                 }
    
    output <- data.frame(    
        sampleInterval_start = sampleInterval_start,
        sampleInterval_end = sampleInterval_end,
        sampleSedColumnDepth_start = fossilSeries$sampleIntervals[,1],
        sampleSedColumnDepth_end = fossilSeries$sampleIntervals[,2],
        sampleMidAge = apply(sampleIntervalAges, 1, mean),
        bioturbInterval_start = fossilSeries$bioturbIntervals[,1], 
        bioturbInterval_end = fossilSeries$bioturbIntervals[,2],
        trueGradientValue = sampleGradientValues,
        isBackgroundSegment = isBackgroundSegment,
        eventID = eventID
        )
    
    return(output)
    }



ageMatchFun <- function(age, eventStartEndTimes){
        # which events does it overlap with
        ageMatch <- (
            age[1] > eventStartEndTimes[,2]  &  
            eventStartEndTimes[,1] > age[2]
            )
        
        # test ageMatch - if none, make NA
        if(any(ageMatch)){
            if(sum(ageMatch) > 1){
                stop(
                    "More than one event matching to a sample! Need to increase intervening background gaps between events"
                    )
                }
            # if one ageMatch, change to ID of which ageMatch
            ageMatch <- which(ageMatch)
        }else{
            ageMatch <- NA
            }
        
        return(ageMatch)
        }