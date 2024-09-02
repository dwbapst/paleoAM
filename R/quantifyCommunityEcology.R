#' Quantify Community Ecology From a Time Series of Simulated Fossil Assemblages
#' 
#' Given the output from functions such as \code{simulateFossilAssemblageSeries}, 
#' this function calculates various ecological measures that might be used as
#' 'recovered' gradient values to compare against the original generating gradient values.

#' @details
#' Quantifies 

#' @inheritParams getSampleDCA The 

#' @param fossilSeries The


#' @return
#' A list containing the input data and arguments, as well as the 
#' DCA score values found  

# @aliases

# @seealso

# @references

# @examples


#' @name quantifyCommunityEcology
#' @rdname quantifyCommunityEcology
#' @export
quantifyCommunityEcology <- function(
        origAbundData, 
        fossilSeries, 
        useTransformedRelAbundance = TRUE,
        projectIntoOrigDCA = TRUE,
        powerRootTransform = 1
    
        #inclusiveDCA = FALSE, 
        #singularDCA = TRUE, 
        #rawDCA = FALSE,
    
        ){
    
    #require(vegan)
  
    ecologyOutList <- list(
        simAbundanceTable = fossilSeries$abundanceTable,
        origAbundanceTable = origAbundData,
        useTransformedRelAbundance = useTransformedRelAbundance,
        projectIntoOrigDCA = projectIntoOrigDCA,
        powerRootTransform = powerRootTransform
        #artificialAbundanceTable = abundanceTable,
        #braycurtisDistMat = bcdist,
        #dcaOut = dcaOut
        )    
    
    # 'singular' DCA is now the only real option
    # if(singularDCA){
    
    scoreDCA1_recovered <- numeric()
    for(i in 1:nrow(fossilSeries$abundanceTable)){
        # use getSampleDCA
        scoreDCA1_recovered[i] <- getSampleDCA(
            simSample = fossilSeries$abundanceTable[i,],
            origAbundData = origAbundData,
            useTransformedRelAbundance = useTransformedRelAbundance,
            projectIntoOrigDCA = projectIntoOrigDCA,
            whichAxes = 1,
            powerRootTransform = powerRootTransform,
            returnDCAforOrigAndSim = FALSE
            )
        }
    
    ecologyOutList$scoreDCA1_recovered <- scoreDCA1_recovered
        
    # }
    
    ######################
    # deprecated ways of getting the DCA value for     
    
    #if(rawDCA){
    #    abundanceTable <- fossilSeries$abundanceTable                        # 1
    #    # Transform the abundance table to relative abundances, 
    #      # and get the pair-wise Bray-Curtis distances.
    #    # recalculate relative abundance table
    #    relAbundanceTable <- t(apply(abundanceTable, 1, function(x) x/sum(x)))
    #    # get bray-curtis distances
    #    bcdist <- vegan::vegdist(relAbundanceTable, method = "bray")
    #    ## Doing a DCA on the Simulated Data
    #    dcaOut <- vegan::decorana(relAbundanceTable)
    #    scoreDCA1_raw <- vegan::scores(dcaOut)[ , 1]
    #    ecologyOutList$scoreDCA1_raw <- scoreDCA1_raw
    #    }    
    
    #if(inclusiveDCA){
    #    # 05-15-21: Need to combine original abundance data 
    #      # into simulated to properly scale DCA (maybe?)
    #      # combine abundance table with original abundance table
    #      # add TEN copies of the original data
    #    abundanceTable <- rbind(fossilSeries$abundanceTable, origAbundData)  # 1
    #    abundanceTable <- rbind(abundanceTable, origAbundData)             # 2
    #    abundanceTable <- rbind(abundanceTable, origAbundData)             # 3
    #    abundanceTable <- rbind(abundanceTable, origAbundData)             # 4
    #    abundanceTable <- rbind(abundanceTable, origAbundData)             # 5
    #    abundanceTable <- rbind(abundanceTable, origAbundData)             # 6
    #    abundanceTable <- rbind(abundanceTable, origAbundData)             # 7
    #    abundanceTable <- rbind(abundanceTable, origAbundData)             # 8
    #    abundanceTable <- rbind(abundanceTable, origAbundData)             # 9
    #    abundanceTable <- rbind(abundanceTable, origAbundData)             # 10
    #    # Transform the abundance table to relative abundances, 
    #      # and get the pair-wise Bray-Curtis distances.
    #    # recalculate relative abundance table
    #    relAbundanceTable <- t(apply(abundanceTable, 1, function(x) x/sum(x)))
    #    # get bray-curtis distances
    #    bcdist <- vegan::vegdist(relAbundanceTable, method = "bray")
    #    ## Doing a DCA on the Simulated Data
    #    dcaOut <- vegan::decorana(relAbundanceTable)
    #    scoreDCA1_inclusive <- vegan::scores(dcaOut)[1:nrow(fossilSeries$abundanceTable), 1]
    #    ecologyOutList$scoreDCA1_inclusive <- scoreDCA1_inclusive
    #    }
        
    return(ecologyOutList)
    }

