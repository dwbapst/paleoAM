#' Summarize Simulation DCA
#' 
#' The

#' @details
#' The

# @inheritParams

#' @param sampleData The

#' @param nEvents The

#' @param backgroundUpperBound The

#' @param eventMagnitudesAbsolute The

#' @param eventLowerBound The

#' @param eventUpperBound The

#' @return
#' The

# @examples


#' @name summarizeSimulationDCA
#' @rdname summarizeSimulationDCA
#' @export
summarizeSimulationDCA <- function(
            sampleData,
            nEvents,
            backgroundUpperBound,
            eventMagnitudesAbsolute,
            eventLowerBound,
            eventUpperBound
            ){

    # summarize a set of simulation data, using their DCA scores 
        # to test if events are significantly different from values 
    # estimated based on background-level simulations (by default, a separately generated sample)

    # find all the unique events that overlap with samples    
    uniqueEvents <- sort(unique(
        sampleData$eventID[!is.na(sampleData$eventID)]
        ))

    # get a list that records all obs DCA values within the interval for each event
    eventDCAvalues <- lapply(uniqueEvents, function(x)
                        sampleData$scoreDCA1_recovered[
                            sampleData$eventID == x & !is.na(sampleData$eventID)
                            ]
                        )

    # measure how many DCA measurements are detected, 
        # ie. outside expected magnitude 
        # for the background during each event
    isEventDetected <- sapply(eventDCAvalues, function(x) 
                            sum(x > backgroundUpperBound)
                            )
        
    # measure how many DCA measurements are detected 
        # within expected magnitude for each event
        # based on simulated 95% envelope
    isDetectedEventWithinExpMagnitude <- sapply(eventDCAvalues, function(x) 
                sum(x > backgroundUpperBound &
                    x > eventLowerBound &
                    x < eventUpperBound
                    )
                )

    # measure how many DCA measurements are detected 
        # within 1 DCA unit for each event
    isDetectedEventWithinOneUnit <- sapply(eventDCAvalues, function(x) 
            sum(x > backgroundUpperBound &
                x > (eventMagnitudesAbsolute - 0.5) &
                x < (eventMagnitudesAbsolute + 0.5)
                )
            )

    # measure how many DCA measurements are detected 
        # within 0.5 DCA unit for each event
    isDetectedEventWithinHalfUnit <- sapply(eventDCAvalues, function(x) 
            sum(x > backgroundUpperBound &
                x > (eventMagnitudesAbsolute - 0.25) &
                x < (eventMagnitudesAbsolute + 0.25)
                )
            )
    
    
    # one more check
    if(any(isDetectedEventWithinExpMagnitude > isEventDetected)){
        stop("More events detected within expected magnitude than TOTAL events just detected!")
        }

    summaryOut <- list()
    summaryOut$backgroundUpperBound <- backgroundUpperBound
    summaryOut$eventLowerBound <- eventLowerBound
    summaryOut$eventUpperBound <- eventUpperBound
    #summaryOut$maxDCAsampled <- maxDCAsampled
    #summaryOut$minDCAsampled <- minDCAsampled
    summaryOut$nEventsDetected <- sum(isEventDetected > 0)
    #summaryOut$nEventsWithinExpectedMagnitude <- sum(isEventWithinExpMagnitude > 0)
    summaryOut$nDetectedEventsWithinExpectedMagnitude <- sum(isDetectedEventWithinExpMagnitude > 0)
    summaryOut$nDetectedEventsWithinOneUnit <- sum(isDetectedEventWithinOneUnit > 0)
    summaryOut$nDetectedEventsWithinHalfUnit <- sum(isDetectedEventWithinHalfUnit > 0)
    summaryOut$eventDCAvaluesList <- eventDCAvalues
    
    return(summaryOut)
    }

