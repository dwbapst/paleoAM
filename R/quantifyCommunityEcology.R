quantifyCommunityEcology <- function(origAbundData, coreRecord, 
    singularDCA = TRUE, inclusiveDCA = FALSE, rawDCA = FALSE){
    
    #require(vegan)
  
    ecologyOutList <- list(
        simAbundanceTable = coreRecord$abundanceTable,
        origAbundanceTable = origAbundData
        #artificialAbundanceTable = abundanceTable,
        #braycurtisDistMat = bcdist,
        #dcaOut = dcaOut
        )    

    if(singularDCA){
        scoreDCA1_singular <- numeric()
        for(i in 1:nrow(coreRecord$abundanceTable)){
            abundanceTable <- rbind(coreRecord$abundanceTable[i,], origAbundData)
            # Transform the abundance table to relative abundances, 
                # and get the pair-wise Bray-Curtis distances.
                # recalculate relative abundance table
            relAbundanceTable <- t(apply(abundanceTable, 1, function(x) x/sum(x)))
            # get bray-curtis distances
            bcdist <- vegan::vegdist(relAbundanceTable, method = "bray")
            ## Doing a DCA on the Simulated Data
            dcaOut <- vegan::decorana(relAbundanceTable)
            scoreDCA1_singular[i] <- vegan::scores(dcaOut)[1, 1]
            }
        
        ecologyOutList$scoreDCA1_singular <- scoreDCA1_singular
        }
    
    if(rawDCA){
        abundanceTable <- coreRecord$abundanceTable                        # 1
        # Transform the abundance table to relative abundances, 
          # and get the pair-wise Bray-Curtis distances.
        # recalculate relative abundance table
        relAbundanceTable <- t(apply(abundanceTable, 1, function(x) x/sum(x)))
        # get bray-curtis distances
        bcdist <- vegan::vegdist(relAbundanceTable, method = "bray")
        ## Doing a DCA on the Simulated Data
        dcaOut <- vegan::decorana(relAbundanceTable)
        scoreDCA1_raw <- vegan::scores(dcaOut)[ , 1]
        ecologyOutList$scoreDCA1_raw <- scoreDCA1_raw
        }    
    
    if(inclusiveDCA){
        # 05-15-21: Need to combine original abundance data into simulated to properly scale DCA (maybe?)
          # combine abundance table with original abundance table
          # add TEN copies of the original data
        abundanceTable <- rbind(coreRecord$abundanceTable, origAbundData)  # 1
        abundanceTable <- rbind(abundanceTable, origAbundData)             # 2
        abundanceTable <- rbind(abundanceTable, origAbundData)             # 3
        abundanceTable <- rbind(abundanceTable, origAbundData)             # 4
        abundanceTable <- rbind(abundanceTable, origAbundData)             # 5
        abundanceTable <- rbind(abundanceTable, origAbundData)             # 6
        abundanceTable <- rbind(abundanceTable, origAbundData)             # 7
        abundanceTable <- rbind(abundanceTable, origAbundData)             # 8
        abundanceTable <- rbind(abundanceTable, origAbundData)             # 9
        abundanceTable <- rbind(abundanceTable, origAbundData)             # 10
        # Transform the abundance table to relative abundances, 
          # and get the pair-wise Bray-Curtis distances.
        # recalculate relative abundance table
        relAbundanceTable <- t(apply(abundanceTable, 1, function(x) x/sum(x)))
        # get bray-curtis distances
        bcdist <- vegan::vegdist(relAbundanceTable, method = "bray")
        ## Doing a DCA on the Simulated Data
        dcaOut <- vegan::decorana(relAbundanceTable)
        scoreDCA1_inclusive <- vegan::scores(dcaOut)[1:nrow(coreRecord$abundanceTable), 1]
        ecologyOutList$scoreDCA1_inclusive <- scoreDCA1_inclusive
        }
        
    return(ecologyOutList)
    }

