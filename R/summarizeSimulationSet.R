#' Summarize a Simulation Set
#' 
#' The

#' @details
#' The

# @inheritParams

#' @param simSummaryTable The

#' @param eventGradientChange_paramValues The

#' @param eventSampleWidth_paramValues The

#' @return
#' The

# @aliases

# @seealso

# @references

# @examples


#' @name summarizeSimulationSet
#' @rdname summarizeSimulationSet
#' @export
summarizeSimulationSet <- function(
            simSummaryTable, 
            eventGradientChange_paramValues, 
            eventSampleWidth_paramValues
            ){

    # summarize a set of simulations, counting how many events are detected
        # and how many fall within an expected magnitude
    # output as a series of tables with event/sample width ratios included
    
    ############################################################    
    propEventsDetected <- (
        simSummaryTable$nEventsDetected / simSummaryTable$nEvents
        )

    #propEventsWithinExpMag <- (
    #    simSummaryTable$nEventsWithinExpectedMagnitude / simSummaryTable$nEvents
    #    )

    propEventsDetectedAndWithinExpMag <- (
        simSummaryTable$nDetectedEventsWithinExpectedMagnitude / simSummaryTable$nEvents
        )
    
    propEventsDetectedAndWithinOneUnit <- (
        simSummaryTable$nDetectedEventsWithinOneUnit / simSummaryTable$nEvents
        )
    
    propEventsDetectedAndWithinHalfUnit <- (
        simSummaryTable$nDetectedEventsWithinHalfUnit / simSummaryTable$nEvents
        )    

    #propDetectedEventsWithinExpMag <- (
    #    simSummaryTable$nDetectedEventsWithinExpectedMagnitude / simSummaryTable$nEventsDetected
    #    )
    # set those with 0 detected events to 0
    #propDetectedEventsWithinExpMag[simSummaryTable$nEventsDetected < 1] <- 0
    
    ###############
    # Checks
    
    if(any(is.na(propEventsDetected))){
        stop("some values of propEventsDetected are NA")
        }
    
    #if(any(is.na(propEventsWithinExpMag))){
    #   stop("some values of propEventsWithinExpMag are NA")
    #   }
    
    if(any(is.na(propEventsDetectedAndWithinExpMag))){
        stop("some values of propEventsDetectedAndWithinExpMag are NA")
        }
    
    if(any(propEventsDetected > 1)){
        stop("some values of proportion of events detected greater than 1??")
        }

    #if(any(propEventsWithinExpMag > 1)){
    #   stop("some values of proportion of events within SD of expected magnitude greater than 1??")
    #   }

    if(any(propEventsDetectedAndWithinExpMag > 1)){
        stop("some values of proportion of detected events within SD of expected magnitude greater than 1??")
        }

    if(any(propEventsDetected < 0)){
        stop("some values of proportion of events detected less than 0??")
        }

    #if(any(propEventsWithinExpMag < 0)){
    #   stop("some values of proportion of events within SD of expected magnitude less than 0??")
    #   }

    if(any(propEventsDetectedAndWithinExpMag < 0)){
        stop("some values of proportion of detected events within SD of expected magnitude less than 0??")
        }
        
    ####################################################
    
    makeSimSetTableTemplate <- function(fill_variable){
        makeSimSetTable(
            fill_variable = fill_variable,
            
            y_variable_from_table = simSummaryTable$eventGradientChange_param,
            y_variable = eventGradientChange_paramValues,
            y_variable_label = "eventGradientChangeScale_",
            
            x_variable_from_table = simSummaryTable$eventSampleWidth_param,
            x_variable = eventSampleWidth_paramValues,
            x_variable_label = "eventSampleWidthRatio_"
            )
        }
        
    #########################################################
            
    propEventsDetected <- makeSimSetTableTemplate(
        fill_variable = propEventsDetected)

    propEventsDetectedAndWithinExpMag <- makeSimSetTableTemplate(
        fill_variable = propEventsDetectedAndWithinExpMag)
    
    propEventsDetectedAndWithinOneUnit <-makeSimSetTableTemplate(
        fill_variable = propEventsDetectedAndWithinOneUnit)
    
    propEventsDetectedAndWithinHalfUnit <-makeSimSetTableTemplate(
        fill_variable = propEventsDetectedAndWithinHalfUnit)    
    
    backgroundUpperBound <- makeSimSetTableTemplate(
        fill_variable = simSummaryTable$backgroundUpperBound)

    eventLowerBound <- makeSimSetTableTemplate(
        fill_variable = simSummaryTable$eventLowerBound)

    eventUpperBound <- makeSimSetTableTemplate(
        fill_variable = simSummaryTable$eventUpperBound)
        
    medianMaxDCA <- makeSimSetTableTemplate(
        fill_variable = simSummaryTable$medianMaxDCA_acrossEvents)
    
    medianPooledDCA <- makeSimSetTableTemplate(
        fill_variable = simSummaryTable$medianPooledDCA)

    ###########################################################
    
    # get the maximum of the medians calculated on pooled DCA values
    maxMedianPooledDCA <- max(simSummaryTable$medianPooledDCA)
    
    simSetSummaryList <- list(
        simSummaryTable = simSummaryTable,
        propEventsDetected = propEventsDetected,
        #propEventsWithinExpMag = propEventsWithinExpMag,
        #propDetectedEventsWithinExpMag = propDetectedEventsWithinExpMag,
        propEventsDetectedAndWithinExpMag = propEventsDetectedAndWithinExpMag,
        propEventsDetectedAndWithinOneUnit = propEventsDetectedAndWithinOneUnit,
        propEventsDetectedAndWithinHalfUnit = propEventsDetectedAndWithinHalfUnit,        
        medianMaxDCA = medianMaxDCA,
        medianPooledDCA = medianPooledDCA,
        maxMedianPooledDCA = maxMedianPooledDCA,
        backgroundUpperBound = backgroundUpperBound,
        eventLowerBound = eventLowerBound,
        eventUpperBound = eventUpperBound
        )
    
    return(simSetSummaryList)
    }

# make the sim set table
makeSimSetTable <- function(
                fill_variable,
                x_variable_from_table,
                y_variable_from_table,
                x_variable,
                y_variable,
                x_variable_label = NULL,
                y_variable_label = NULL
                ){
        
    ##################################
    tableOut <- matrix(NA, 
        length(x_variable), 
        length(y_variable)
        )
    
    if(!is.null(x_variable_label)){
        rownames(tableOut) <- paste0(x_variable_label, x_variable)
        }

    if(!is.null(y_variable_label)){
        colnames(tableOut) <- paste0(y_variable_label, y_variable)        
        }

    matchCoordinates <- cbind(
        match(x_variable_from_table, x_variable),
        match(y_variable_from_table, y_variable)
        )
    
    if(length(dim(tableOut)) != 2){
        stop("Table isn't two-dimensional")
        }
    
    if(any(is.na(matchCoordinates))){
        stop("NA produced in coordinate matching")
        }
    
    for(i in 1:nrow(matchCoordinates)){
        xCoord <- matchCoordinates[i,1]
        yCoord <- matchCoordinates[i,2]
        filler <- fill_variable[i]
        if(is.list(filler)){
            filler <- unlist(filler)
            }
        tableOut[xCoord, yCoord] <- filler
        }
    
    return(tableOut)
    }


# add to the simSummaryTable
addToSimSumTable <- function(simSummaryTable, simSummaryOut){

    # assign output to values in simSummaryTable
    simSummaryTable$backgroundUpperBound[i] <- simSummaryOut$backgroundUpperBound
    simSummaryTable$eventLowerBound[i] <- simSummaryOut$eventLowerBound
    simSummaryTable$eventUpperBound[i] <- simSummaryOut$eventUpperBound

    simSummaryTable$nEventsDetected[i] <- simSummaryOut$nEventsDetected
    #simSummaryTable$nEventsWithinExpectedMagnitude[i] <- simSummaryOut$nEventsWithinExpectedMagnitude
    simSummaryTable$nDetectedEventsWithinExpectedMagnitude[i] <- simSummaryOut$nDetectedEventsWithinExpectedMagnitude
    simSummaryTable$nDetectedEventsWithinOneUnit[i] <- simSummaryOut$nDetectedEventsWithinOneUnit 
    simSummaryTable$nDetectedEventsWithinHalfUnit[i] <- simSummaryOut$nDetectedEventsWithinHalfUnit     

    return(simSummaryTable)
    }


# add data to simSetList
addToSimSetList <- function(simSetList, simSummaryOut){

    # save list of DCA values from each event to grand-list for whole set
    simSetList$simSetDCAvalues[[i]] <- simSummaryOut$eventDCAvaluesList

    simSetList$medianDCA[[i]] <- sapply(simSetList$simSetDCAvalues[[i]] , function(x) 
        median(x, na.rm = TRUE)
       )    
    simSetList$meanDCA[[i]] <- sapply(simSetList$simSetDCAvalues[[i]] , function(x) 
        mean(x, na.rm = TRUE)
        )
    simSetList$maxDCA[[i]] <- sapply(simSetList$simSetDCAvalues[[i]] , function(x) 
        max(x, na.rm = TRUE)
        )
    
    return(simSetList)
    }
