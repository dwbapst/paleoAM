#' Plot the Recovered DCA Values of a Fossil Assemblage Series 
#' 
#' Makes a plot of the simulated generating gradient over time along with the
#' recovered gradient values in the same plot for a given simulated time series 
#' of fossil assemblages, particularly those from
#' \code{\link{simulateFossilAssemblageSeries}}.

#' @details
#' The function will generally only be run on the output 
#' from \code{\link{simulateFossilAssemblageSeries}}, although the function is 
#' written so that the necessary elements can be provided separately.

# @inheritParams 

#' @param ... This function takes either the output from 
#' \code{\link{simulateFossilAssemblageSeries}} or requires specifying 
#' the three arguments \code{simTimeVar}, a data-frame of time-steps, 
#' sedimentary width and gradient values for a time-series simulation); 
#' \code{gradientRecovered}, the recovered gradient values; and
#' \code{sampleAge}, the age of individual samples.
 
#' @param colSimGenerating What color should be used for the generating 
#' ("true") gradient values?

#' @param colSimRecovered What color should be used for the recovered 
#' gradient values calculated from the simulated data?

#' @return
#' Returns nothing at all. Just a plot. That's all!

# @aliases

# @seealso

# @references

# @examples


#' @name plotFossilAssemblageSeriesDCA
#' @rdname plotFossilAssemblageSeriesDCA
#' @export
plotFossilAssemblageSeriesDCA <- function(
            ..., 
            colSimGenerating = "black", 
            colSimRecovered = "navy"
            ){

    arguments <- list(...)
    #

    if(length(arguments) == 1){

        simFossilAssemblageSeriesOut <- arguments[[1]]
        plotFossilAssemblageSeriesDCA_sepVariables(
            simTimeVar = simFossilAssemblageSeriesOut$simTimeVar, 
            gradientRecovered = simFossilAssemblageSeriesOut$sampleProperties$scoreDCA1_recovered,
            sampleAge = simFossilAssemblageSeriesOut$sampleProperties$sampleMidAge,
            colSimGenerating = colSimGenerating, 
            colSimRecovered = colSimRecovered
            )

    }else{
    
        if(length(arguments) == 3){
            
            plotFossilAssemblageSeriesDCA_sepVariables(
                simTimeVar = arguments$simTimeVar, 
                gradientRecovered = arguments$gradientRecovered, 
                sampleAge = arguments$sampleAge,
                colSimGenerating = colSimGenerating, 
                colSimRecovered = colSimRecovered
                )
            
        }else{
            stop(paste0(
                "Wrong number of input arguments for plotting: use either output from\n",
                "simulateFossilAssemblageSeries or use the\n three elements 'simTimeVar', 'gradientRecovered', 'sampleAge'"
                ))
            }
        }
    }

plotFossilAssemblageSeriesDCA_sepVariables <- function(
    simTimeVar, 
    gradientRecovered, 
    sampleAge,
    colSimGenerating = "black", 
    colSimRecovered = "navy"
    ){
    
    # calculate y limits
    DCA1lump <- c(simTimeVar$gradientValue, gradientRecovered)
    DCA1_ylims <- c(min(DCA1lump),max(DCA1lump))
    DCA1_ylims[2] <- DCA1_ylims[2] + ((DCA1_ylims[2] - DCA1_ylims[1])/3)
    
    #par(mar = c(5,5,6,1))
    
    plot(simTimeVar$timestep, simTimeVar$gradientValue,
         main = "Original vs Measured Values for DCA-1 Score",
         #main = paramTitle,
         xlab = "Age (time-steps)",
         ylab = "DCA Score 1",
         lty = 1, lwd = 1,
         ylim = DCA1_ylims,
         col = colSimGenerating,
         type = "l"
         )
    
    graphics::points(sampleAge, gradientRecovered, 
           cex = 1, pch = 16, 
           col = colSimRecovered
           )
    
    graphics::legend(x = "top",
           legend = c("Original", "Measured"),
           lty = c(1,3), 
           lwd = c(1,3), 
           ncol = 2,
           col = c(colSimGenerating, colSimRecovered), 
           bg = "white"
           )
    
    }
