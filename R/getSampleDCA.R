getSampleDCA <- function(
            simPickedSample,
            origAbundData,
            useRelAbundance = TRUE,
            returnDCAforOrigAndSim = FALSE,
            whichAxes = 1
            ){
        # take a simulated sample, combine with original data,
            # apply DCA and get the DCA-1 value
        
        # combined picked sample with original abundance data
        abundanceTable <- rbind(origAbundData, simPickedSample)

        if(useRelAbundance){        
            # Transform the abundance table to relative abundances, 
                # and get the pair-wise Bray-Curtis distances.
                # recalculate relative abundance table
            relAbundanceTable <- t(apply(abundanceTable, 1, function(x) x/sum(x)))
            ## Doing a DCA on the Simulated Data
            dcaOut <- vegan::decorana(relAbundanceTable)
        }else{
            ## Doing a DCA on the Simulated Data
            dcaOut <- vegan::decorana(abundanceTable)
            }
        
        if(returnDCAforOrigAndSim){
            # get DCA value for sim assemblages AND empirical data
            output <- vegan::scores(dcaOut)[ , whichAxes]
        }else{
            # get DCA value for the simulated mixed assemblage
            output <- vegan::scores(dcaOut)[nrow(abundanceTable), whichAxes]
            }
        
        # get bray-curtis distances
        #bcdist <- vegan::vegdist(relAbundanceTable, method = "bray")

        return(output)
        }
