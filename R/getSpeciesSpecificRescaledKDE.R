#' This is a function for Fitting a KDE to a specific species in Community Ecology Data

#' @details
#' This function fits a KDE to the abundance data of a particular species from community data given some environmental gradient.

#' @param gradientOrigDCA The environmental gradient along which abundance varies, which you are fitting a KDE to.
 
#' @param origAbundData The abundance data of the data you wish to model the abundance of.
 
#' @param abundanceFloorRatio The minimum value for the abundance in a given interval along the gradient -- a probably arbitration value that is set to 0.5 by default.
 
#' @param nBreaksGradientHist The default is 20. Twenty what they asked? Twenty something.

#' @param modeledSiteAbundance The number of abundances the relative abundances will by multiplied by to formulate the KDE. The default is 10000.

# @param reportTests Should the result of tests checking the KDE results be reported to the terminal? Default is /code{FALSE}.

#' @return
#' A list containing the KDEs describing change in abundance for each species across the specified gradient.

#' @seealso
#' \code{\link{getProbOccViaPresAbs}}

#' @references

#' @examples
#' 
#' # load data
#' data(gulfOfAlaska_ShannonEtAl)
#' 
#' gradientOrigDCA <- getSpeciesSpecificRescaledKDE(
#'     gradientOrigDCA = DCA1_GOA, 
#'     origAbundData = abundData_GOA, 
#'     abundanceFloorRatio = 0.5, 
#'     nBreaksGradientHist = 20, 
#'     modeledSiteAbundance = 10000
#'     )
#'     


probSpeciesOccur <- getProbOccViaPresAbs(
    gradientOrigDCA = DCA1_GOA, 
    origAbundData = abundData_GOA
    )

#' @name getSpeciesSpecificRescaledKDE
#' @rdname getSpeciesSpecificRescaledKDE
#' @export
getSpeciesSpecificRescaledKDE <- function(
        gradientOrigDCA, 
        origAbundData,
        abundanceFloorRatio = 0.5,
        nBreaksGradientHist = 20,
        modeledSiteAbundance = 10000
        ){
    
    # Species-Specific Kernel Density Estimates on the Original Empirical DCA-1 Gradient
    # Obtain species KDEs for simulating from an empirical DCA-1 gradient
    
    if(length(gradientOrigDCA) != nrow(origAbundData)){
        stop("DCA and abundance data seem to not agree on number of samples/sites")
        }
    
    gradientHist <- hist(
        gradientOrigDCA, 
        breaks = nBreaksGradientHist, 
        plot = FALSE
        #,main = "Number of Samples in each Binned Segment of Gradient"
        #,xlab = "Gradient (DCA-1 Values)"
        )
    
    histBreaks <- gradientHist$breaks
    histSampleCounts <- gradientHist$counts
    
    if(any(histSampleCounts<1)){
        stop(
            "Cannot have regions of gradient with specified number of breaks with no samples. Try fewer breaks."
            )
        }
    
    # Convert absolute abundance data to relative abundance, 
    # multiply by arbitrarily large number (modeledSiteAbundance)
        # to make the abundances into whole numbers
    
    # convert to relative abundances, then multiply each by modeledSiteAbundance
    origAbundData_RelScale <- t(apply(origAbundData, 1, function(x) x/sum(x)))
    origAbundData_RelScale <- as.data.frame(round(
        origAbundData_RelScale * modeledSiteAbundance
        ))
    
    # get the smallest non-zero abundance
    smallestRelAbund <- unique(as.vector(unlist(origAbundData_RelScale)))
    smallestRelAbund <- min(smallestRelAbund[smallestRelAbund > 0])
    
    # Add X to all sites so all species have >0 for all sites
        # so each species appears once in every sample
    # add smallestRelAbund * abundanceFloorRatio
        # default for abundanceFloorRation is 0.5
        # so default is every species treated as 
            # occurring at half of the smallest relative abundance 
            # reported across the entire data set
    abundFloorModifier <- ceiling(abundanceFloorRatio * smallestRelAbund)
    
    origAbundData_RelScale <- origAbundData_RelScale + abundFloorModifier
    
    # Take each gradient * abundance # for each species. 
    # Replicate gradient values for *each* specimen 
    
    gradientRepAbund <- apply(
        origAbundData_RelScale, 2, 
          function(y){
            #if(any(y<1)){
            #    stop(paste0("This vector has a zero: ", y))
            #    }
            #if(any(is.na(y))){
            #    stop(paste0("This vector has an NA: ", y))
            #   }
            combMat <- cbind(gradientOrigDCA, y)
            #print(combMat)
            unlist(apply(combMat, 1, 
                function(x) rep(x = x[1], times = x[2]))
                )
            }
        )
    
    # count abundance in each bin from the hist breaks
    abundanceHistBreaks <- lapply(gradientRepAbund, function(x) 
        hist(x, breaks = histBreaks, plot = FALSE
            )$counts
        )
    
    # These repeated gradient values don't account 
    # for the non-random sampling of the gradient itself. 
    # We need to bin, and normalize abundances relative to the sampling of the gradient
    
    # 05-17-22
    # need to condition abundances on the number of samples **they were sampled in**
    
    # convert abundance data to presence data
    origPresenceData <- origAbundData
    origPresenceData[origPresenceData > 0] <- 1
    # Count the number of presences across the gradient
        # convert to DCA values for applying hist() to get
        # frequency within bins
    gradientRepPresence <- apply(
        origPresenceData, 2,
        function(x) gradientOrigDCA[x > 0]                         
        )
    # count sites each species occurs at in each bin
    gradientPresenceHistBreaks <- lapply(gradientRepPresence, function(x) 
        hist(x, 
            breaks = histBreaks, 
            plot = FALSE
            )$counts + 1
        )    

    # Now get ratio of gradientRepAbund to histSampleCounts 
    # This makes the number of occurrences proportional to sampling.
        # 05-17-22
        # needs to be made proportional to the number of samples that
        # each individual species is actually observed in, not all samples
    gradientSampCorAbundHistBreaks <- lapply(1:length(abundanceHistBreaks),
        function(x) ceiling(abundanceHistBreaks[[x]]/gradientPresenceHistBreaks[[x]])
        )
    names(gradientSampCorAbundHistBreaks) <- names(abundanceHistBreaks)
 
    # now remake gradientRepAbund
    gradientRepSampAbund <- lapply(gradientSampCorAbundHistBreaks, 
       function(y) unlist(
           apply(cbind(gradientHist$mids,y), 1, 
                 function(x) rep(x = x[1], times = x[2]))
           )
        )
    
    #hist(gradientRepAbund[[1]])
    # apply(abundData, 2, max)
    # plot(gradient, abundData$EpistominellaExigua)
    # hist(gradientRepAbund$EpistominellaExigua)
    # hist(gradientRepAbund$BuliminellaTenuata)
    
    # Now fit the KDE to each species' abundance curve
    kdePerTaxa <- lapply(gradientRepSampAbund, density, 
        bw = diff(gradientHist$mids)[1] )

    # Get maximum and min abundance for each species 
    # from the sampling normalized gradientSamplingCorrAbund 
    
    # Get mean relative abundance as well. 
    # Assign max/min and mean relative abundance to
        # elements stored within the list object from
        # the kernel density estimation function

    for(i in 1:length(kdePerTaxa)){
        
        #kdePerTaxa[[i]]$meanRelAbund <- mean(origAbundData_RelScale[,i])
        sampCorAbund_taxon <- gradientSampCorAbundHistBreaks[[i]]
        kdePerTaxa[[i]]$maxSampCorrAbund <- max(sampCorAbund_taxon)
        kdePerTaxa[[i]]$minSampCorrAbund <- min(sampCorAbund_taxon)
        
        if(min(sampCorAbund_taxon) < 0){
            stop("there seems to be bad minimum abundance measures")
            }
        if(max(sampCorAbund_taxon) < 0){
            stop("there seems to be bad maximum abundance measures")
            }
        }
    
    # Then we scale area under each KDE to the maximum observed abundance of each KDE.
        # Could we use something else like median abundance?
        # No, because then taxa are overly abundant in a small number of samples 
    # Alternatives are definitely something to ponder for future work.

    rescaleFunctKDE <- function(kdeEst, 
        modeledSiteAbundance){
        
        # min is rescaled to zero 
            # This may be problematic if some species have
            # positive abundance across the whole gradient
        # scale height to the max relative abundance
            # instead of scaling the height of the KDE (max) to the 
            # sampling corrected expected abundance calculated the KDE
        
        # scale to min of 0
        kdeEst$y <- kdeEst$y - min(kdeEst$y)
        # scale to max of 1
        kdeEst$y <- kdeEst$y / max(kdeEst$y)
        
        # scale density height to 0 to max scaled abundance + 1
        
        #if(reportTests){
        #    message( paste0("kdeEst$max is ", kdeEst$max) )
        #    }

        kdeEst$y <- kdeEst$y * (
            kdeEst$maxSampCorrAbund - kdeEst$minSampCorrAbund
            ) 
        
        # shift the bottom of the KDE upwards 
            # so at min sampling-corrected abundance (which might not be zero)
        
        kdeEst$y <- kdeEst$y + kdeEst$minSampCorrAbund
        
        # check to make sure not present
        if(any(kdeEst$y < 0)){
            stop("Negative abundances in KDE??")
            }
        
        # scale density to the modeled site abundance
            # to give expected relative abundance
        kdeEst$y <- kdeEst$y / modeledSiteAbundance
        
        return(kdeEst)
        }
    
    kdeRescaled <- lapply(
        kdePerTaxa, 
        rescaleFunctKDE, 
        modeledSiteAbundance = modeledSiteAbundance
        )
    
    return(kdeRescaled)
    }

