#' Measure the Duration of a Transition Period from Recovered Sequence of Fossil Assemblages with DCA-1 Scores
#' 
#'  

#' @details


#' @param simRecord 
 
#' @param bgUpperEnvelope 

#' @param eventLowerEnvelope 

#' @param returnAsAge 

#' @param trueEventDuration 

#' @param plot 

#' @return

#' @seealso

#' @references

#' @examples

#' @name getRecoveredTransitionDuration
# @rdname
# @aliases
#' @export
getRecoveredTransitionDuration <- function(
        simRecord, bgUpperEnvelope, 
        eventLowerEnvelope = NULL,
        returnAsAge = FALSE, 
        trueEventDuration = NA, 
        plot = FALSE
        ){
    # find when the record deviates from background
            
    if(all(simRecord$sampleProperties$scoreDCA1_singular < bgUpperEnvelope)){
        stop("No value exceeds bgUpperEnvelope")
        }
    if(all(simRecord$sampleProperties$scoreDCA1_singular > bgUpperEnvelope)){
        stop("No value is below bgUpperEnvelope")
        }        
        
    startTran <- which(simRecord$sampleProperties$scoreDCA1_singular > bgUpperEnvelope)[1]
    
    if(is.null(eventLowerEnvelope)){
        endTran <- which(simRecord$sampleProperties$scoreDCA1_singular ==
                max(simRecord$sampleProperties$scoreDCA1_singular))[1]
    }else{
        if(all(simRecord$sampleProperties$scoreDCA1_singular < eventLowerEnvelope)){
            stop("No value exceeds eventLowerEnvelope")
            }
        
        endTran <- which(simRecord$sampleProperties$scoreDCA1_singular >
                eventLowerEnvelope)[1]
        }
    
    if(is.na(startTran)){stop("Could not find startTran")}
    if(is.na(endTran)){stop("Could not find endTran")}
    
    if(startTran > endTran){
        stop("Oh no the transition ends before it starts!")
        }
    
    if(startTran == endTran){
        tranDurationRatio <- 0
    }else{
        tranDuration <- simRecord$sampleProperties$sampleMidAge[c(startTran, endTran)] 
        tranDuration <- abs(diff(tranDuration))
        tranDurationRatio <- tranDuration / simRecord$implicitParameters$eventDuration
        }
    if(returnAsAge){
        tranDurationOut <- tranDurationRatio * trueEventDuration
    }else{
        # return as ratio
        tranDurationOut <- tranDurationRatio
        }

    if(plot){
        plot(type = "b", 
            x = simResult$sampleProperties$sampleMidAge, 
            y=simResult$sampleProperties$scoreDCA1_singular,
            xlab = "Time (time-steps)",
            ylab = "DCA-1 Gradient",
            main = "")
        abline(v=sampleAge[c(startTran, endTran)] )
        text(x = 0, y = 2, pos = 4, paste0("Tran. Dur. =",
            round(tranDurationOut, 2))
            )
        }    
    
    return(tranDurationOut)
    }

