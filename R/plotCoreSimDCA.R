
plotCoreSimDCA <- function(
                    simTimeVar, scoreDCA1, sampleAge,
                    colSimTrue = "black", 
                    colSimCalc = "navy"
                    ){

    # calculate y limits
    DCA1lump <- c(simTimeVar$gradientValue, scoreDCA1)
    DCA1_ylims <- c(min(DCA1lump),max(DCA1lump))
    DCA1_ylims[2] <- DCA1_ylims[2] + ((DCA1_ylims[2] - DCA1_ylims[1])/3)
    
    #par(mar = c(5,5,6,1))
    
    plot(simTimeVar$year, simTimeVar$gradientValue,
         main = "Original vs Measured Values for DCA-1 Score",
         #main = paramTitle,
         xlab = "Age (years)",
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
