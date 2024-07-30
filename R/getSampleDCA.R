#' @details

#' @inheritParams

#' @param

#' @return

#' @aliases

#' @seealso

#' @references

#' @examples

#' @param simPickedSample 
#'
#' @param origAbundData 
#' @param useTransformedRelAbundance 
#' @param projectIntoOrigDCA 
#' @param returnDCAforOrigAndSim 
#' @param whichAxes 
#' @param powerRootTransform 
#'
#' @name
#' @rdname
#' @export


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
                abundanceTable <- powerRootTransformFun(abundanceTable, 
                    powerRootTransform = powerRootTransform)   
                }
            dcaOut <- vegan::decorana(abundanceTable)
            
            newSample <- rbind(origAbundData,simPickedSample)
            if(useTransformedRelAbundance){        
                newSample <- powerRootTransformFun(newSample, 
                    powerRootTransform = powerRootTransform)   
                }

            newSampleDCA <- vegan:::predict.decorana(
                object = dcaOut,
                newdata = newSample,
                rank = 1, type = "sites"
                )
            
            if(returnDCAforOrigAndSim){
                # get DCA value for sim assemblages AND empirical data
                output <- newSampleDCA[ , 1]
            }else{
                # get DCA value for the simulated mixed assemblage
                output <- newSampleDCA[nrow(newSample), 1]
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
                abundanceTable <- powerRootTransformFun(abundanceTable, 
                    powerRootTransform = powerRootTransform)  
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

powerRootTransformFun <- function(abundTable, powerRootTransform){
    # convert to relative abundance
    out <- apply(abundTable, 1, function(x) x/sum(x))
    out <- t(out) ^ (1/powerRootTransform)  
    return(out)
    }

