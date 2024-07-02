getSampleDCA <- function(
            simPickedSample,
            origAbundData,
            useTransformedRelAbundance = TRUE,
            projectIntoOrigDCA = TRUE,
            returnDCAforOrigAndSim = FALSE,
            whichAxes = 1,
            powerRootTransform = 1
            ){

        if(projectIntoOrigDCA){
                
            if(length(whichAxes) != 1){
                stop("Cannot return any axis other than DCA-1 when 'projectIntoOrigDCA = TRUE' ")
                }
            if(whichAxes != 1){
                stop("Cannot return any axis other than DCA-1 when 'projectIntoOrigDCA = TRUE' ")
                }
            
            abundanceTable <- origAbundData
            ## Doing a DCA just on the original data
            if(useTransformedRelAbundance){
                # root transform
                abundanceTable <- t(apply(abundanceTable, 1, function(x) 
                    x/sum(x)))^(1/powerRootTransform)  
                }
            dcaOut <- vegan::decorana(abundanceTable)
            
            newSample <- rbind(origAbundData,simPickedSample)
            if(useTransformedRelAbundance){        
                newSample <- t(apply(newSample, 1, function(x) 
                    x/sum(x)))^(1/powerRootTransform)  
                }
            newSample <- vegan:::predict.decorana(
                object = dcaOut,
                newdata = newSample,
                rank = 1, type = "sites"
                )
            
            if(returnDCAforOrigAndSim){
                # get DCA value for sim assemblages AND empirical data
                output <- newSample[ , 1]
            }else{
                # get DCA value for the simulated mixed assemblage
                output <- newSample[nrow(newSample), 1]
                }
            
        }else{
            # take a simulated sample, combine with original data,
                # apply DCA and get the DCA-1 value

            # combined picked sample with original abundance data
            abundanceTable <- rbind(origAbundData, simPickedSample)
            #    
            if(useTransformedRelAbundance){        
                # Transform the abundance table to relative abundances, 
                    # (and get the pair-wise Bray-Curtis distances... wait no)
                    # recalculate relative abundance table
                #relAbundanceTable <- t(apply(abundanceTable, 1, function(x) x/sum(x)))
                # root transform?
                abundanceTable <- t(apply(abundanceTable, 1, function(x) 
                    x/sum(x)))^(1/powerRootTransform)  
                }
            ## Doing a DCA on the Simulated Data
            dcaOut <- vegan::decorana(abundanceTable)
            dcaOut <- vegan::scores(dcaOut)
            
            if(returnDCAforOrigAndSim){
                # get DCA value for sim assemblages AND empirical data
                output <- dcaOut[ , whichAxes]
            }else{
                # get DCA value for the simulated mixed assemblage
                output <- dcaOut[nrow(abundanceTable), whichAxes]
                }
            }
        
        # get bray-curtis distances
        #bcdist <- vegan::vegdist(relAbundanceTable, method = "bray")

        return(output)
        }
