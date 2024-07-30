#' @details

#' @inheritParams

#' @param

#' @return

#' @aliases

#' @seealso

#' @references

#' @examples

#' @param fossilSeries 
#'
#' @param minPickedSampleSize 
#' @param maxPickedSampleSize 
#' @param maxDominance 
#'
#' @name
#' @rdname
#' @export

checkFossilSeriesAbundanceTable <- function(
            fossilSeries,
            minPickedSampleSize,
            maxPickedSampleSize,
            maxDominance
            ){
    
    abundanceTable <- fossilSeries$abundanceTable
    sampleIntervals <- fossilSeries$sampleIntervals
    bioturbIntervals <- fossilSeries$bioturbIntervals
    
    # need to throw out samples with too few individuals
    sufficientSampleSize <- rowSums(abundanceTable) >= minPickedSampleSize
    if(sum(sufficientSampleSize) > 2) {
        if(sum(!sufficientSampleSize)!=0){
            message(paste0(sum(!sufficientSampleSize),
                           " samples dropped due to too few specimens sampled"))
            }
        abundanceTable <- abundanceTable[sufficientSampleSize, ]
        sampleIntervals <- sampleIntervals[sufficientSampleSize, ]
        bioturbIntervals <- bioturbIntervals[sufficientSampleSize, ]
        nSamples <- nrow(abundanceTable)
    }else{
        stop("No samples obtained had a sufficient number of specimens for study")
        }
    
    # down-sample samples with too many individuals (>500)
    largeSamples <- rowSums(abundanceTable) > maxPickedSampleSize
    while(any(largeSamples)){
        # use a splitter (cut them in half)
        for(i in which(largeSamples)){
            bigSample <- abundanceTable[i,]  
            # un-table() the community abundance data
            bigSampleList<- rep(names(bigSample), bigSample) 
            # sample to half size
            halfSize <- round(sum(bigSample)/2)
            splitSample <- sample(x = bigSampleList, replace = FALSE, size = halfSize)
            splitSample <- table(splitSample)
            for(j in 1:length(splitSample)){
                speciesMatch <- which(speciesNames == names(splitSample)[j])
                abundanceTable[i, speciesMatch] <- splitSample[j]
                }
            }
        largeSamples <- rowsum(abundanceTable) > maxPickedSampleSize
        }
    
    # test for over-dominated samples
    # get relative abundance table 
    relAbundanceTable <- t(apply(abundanceTable, 1, function(x) x/sum(x)))
    # test for overly-dominated samples (abundance > maxDominance) 
    dominatedSamples <- apply(relAbundanceTable, 1, function(x) any(x > maxDominance))
    if(sum(!dominatedSamples) > 2) {
        if(sum(dominatedSamples) != 0){
            message(paste0(sum(dominatedSamples), 
                           " samples dropped due to over-dominance of a single taxon (>",
                           maxDominance, "%)"))
            }
        abundanceTable <- abundanceTable[!dominatedSamples, ]
        sampleIntervals <- sampleIntervals[!dominatedSamples, ]
        bioturbIntervals <- bioturbIntervals[!dominatedSamples, ]        
        
    }else{
        stop("All samples obtained were over-dominated (>80%) by a single taxon")
        }
    
    fossilSeries$abundanceTable <- abundanceTable
    fossilSeries$sampleIntervals <- sampleIntervals
    fossilSeries$bioturbIntervals <- bioturbIntervals
    
    return(fossilSeries)
    }
