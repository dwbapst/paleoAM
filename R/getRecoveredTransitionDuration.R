#' Measure the Duration of a Transition Period from Recovered Sequence of Fossil Assemblages with DCA-1 Scores
#' 
#' How long did a transition proceed from background to peak 'event' v

#' @details
#' The envelope values can be calculated different ways, or even picked arbitrarily by the user. For example, \code{bgUpperEnvelope} is the upper envelope on what is considered a background value for a gradient value derived from the assemblage (for example, an ordination score). One way a user could calculate \code{bgUpperEnvelope} would be to repeatedly simulate assemblages at the background value, calculate their apparent gradient value and estimate a 0.95 or 0.975 quantile. This can be done easily with function \code{\link{simulateGradientQuantile}}.

#' @param simRecord A simulated fossil record with assemblage change across multiple time-steps, with sedimentary thickness modeled.
 
#' @param bgUpperEnvelope The upper envelope on what is considered a background value for a gradient value derived from the assemblage.

#' @param eventLowerEnvelope The lower envelope on what is considered an event value for a gradient value derived from the assemblage.

#' @param returnAsAge Should the estimated duration of the transition be returned as a duration in time-units? If \code{FALSE} (the default), the value is instead returned as a ratio relative to the true event duration.

#' @param trueEventDuration The true duration of the event. This must be provided by the user if \code{returnAsAge = TRUE} to calculate the duration of the transition interval in simulation time-units.

#' @param plot Should the data be plotted with the estimated transition interval on it, for visual checking? 

#' @return
#' A single value, reflecting (by default) a ratio of transition duration over the event duration. Can be modified with argument \code{returnAsAge}.

#' @seealso
#' \code{\link{simulateGradientQuantile}}

# @references

# @examples

#' @name getRecoveredTransitionDuration
# @rdname
# @aliases
#' @export
getRecoveredTransitionDuration <- function(
        simRecord, 
        bgUpperEnvelope, 
        eventLowerEnvelope = NULL,
        returnAsAge = FALSE, 
        trueEventDuration = NA, 
        plot = FALSE
        ){
    # find when the record deviates from background
            
    if(all(simRecord$sampleProperties$scoreDCA1_recovered < bgUpperEnvelope)){
        stop("No value exceeds bgUpperEnvelope")
        }
    if(all(simRecord$sampleProperties$scoreDCA1_recovered > bgUpperEnvelope)){
        stop("No value is below bgUpperEnvelope")
        }        
        
    startTran <- which(simRecord$sampleProperties$scoreDCA1_recovered > bgUpperEnvelope)[1]
    
    if(is.null(eventLowerEnvelope)){
        endTran <- which(simRecord$sampleProperties$scoreDCA1_recovered ==
                max(simRecord$sampleProperties$scoreDCA1_recovered))[1]
    }else{
        if(all(simRecord$sampleProperties$scoreDCA1_recovered < eventLowerEnvelope)){
            stop("No value exceeds eventLowerEnvelope")
            }
        
        endTran <- which(simRecord$sampleProperties$scoreDCA1_recovered >
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
            x = simRecord$sampleProperties$sampleMidAge, 
            y=simRecord$sampleProperties$scoreDCA1_recovered,
            xlab = "Time (time-steps)",
            ylab = "DCA-1 Gradient",
            main = "")
        graphics::abline(v=simRecord$sampleProperties$sampleMidAge[
            c(startTran, endTran)] 
            )
        graphics::text(x = 0, y = 2, pos = 4, paste0("Tran. Dur. =",
            round(tranDurationOut, 2))
            )
        }    
    
    return(tranDurationOut)
    }
