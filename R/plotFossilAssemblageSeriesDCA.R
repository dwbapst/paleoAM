plotFossilAssemblageSeriesDCA <- function(
            ..., colSimTrue = "black", colSimCalc = "navy"
            ){
    arguments <- list(...)
    #
    if(length(arguments) == 1){
        simFossilAssemblageSeriesOut <- arguments[[1]]
        plotFossilAssemblageSeriesDCA_sepVariables(
            simTimeVar = simFossilAssemblageSeriesOut$simTimeVar, 
            scoreDCA1 = simFossilAssemblageSeriesOut$ecology$scoreDCA1_singular,
            sampleAge = simFossilAssemblageSeriesOut$sampleProperties$sampleMidAge,
            colSimTrue = colSimTrue, 
            colSimCalc = colSimCalc
            )
    }else{
    if(length(arguments) == 3){
        plotFossilAssemblageSeriesDCA_sepVariables(
            simTimeVar = arguments$simTimeVar, 
            scoreDCA1 = arguments$scoreDCA1, 
            sampleAge = arguments$sampleAge,
            colSimTrue = colSimTrue, 
            colSimCalc = colSimCalc
            )
    }else{
        stop(
            "Wrong number of input arguments for plotting: use either output from\n simulateFossilAssemblageSeries or use the\n three elements 'simTimeVar', 'scoreDCA1', 'sampleAge'"
            )}
        }
    }

plotFossilAssemblageSeriesDCA_sepVariables <- function(
                    simTimeVar, scoreDCA1, sampleAge,
                    colSimTrue = "black", 
                    colSimCalc = "navy"
                    ){

    # calculate y limits
    DCA1lump <- c(simTimeVar$gradientValue, scoreDCA1)
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
         col = colSimTrue,
         type = "l"
         )
    
    points(sampleAge, scoreDCA1, 
           cex = 1, pch = 16, 
           col = colSimCalc
           )
    
    legend(x = "top",
           legend = c("Original", "Measured"),
           lty = c(1,3), 
           lwd = c(1,3), 
           ncol = 2,
           col = c(colSimTrue, colSimCalc), 
           bg = "white"
           )
    
    }
