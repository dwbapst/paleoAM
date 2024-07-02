#' This is a function for Fitting a KDE to a specific species in Community Ecology Data

#' @details
#' This function fits a KDE to the abundance data of a particular species from community data given some environmental gradient.

#' @param
#' gradientOrigDCA The environmental gradient along which abundance varies, which you are fitting a KDE to.

#' origAbundData The abundance data of the data you wish to model the abundance of.

#' abundanceFloorRatio The minimum value for the abundance in a given interval along the gradient -- a probably arbitration value that is set to 0.5 by default.
 
#' nBreaksGradientHist The default is 20. Twenty what they asked? Twenty something.

#' modeledSiteAbundance The number of abundances the relative abundances will by multiplied by to formulate the KDE. The default is 10000.

#' reportTests Should the result of tests checking the KDE results be reported to the terminal? Default is /code{FALSE}.

#' # example
#' 
#' # load data
#' data(gulfOfAlaska_ShannonEtAl)
#' 
#' gradientOrigDCA <- getSpeciesSpecificRescaledKDE(
#'     gradientOrigDCA, 
#'     origAbundData, 
#'     abundanceFloorRatio = 0.5, 
#'     nBreaksGradientHist = 20, 
#'     modeledSiteAbundance = 10000, 
#'     reportTests = FALSE)
#'     

getSpeciesSpecificRescaledKDE <- function(
        gradientOrigDCA, 
        origAbundData,
        abundanceFloorRatio = 0.5,
        nBreaksGradientHist = 20,
        modeledSiteAbundance = 10000,
        reportTests = FALSE
        ){
    
    # Species-Specific Kernal Density Estimates on the Original Empirical DCA-1 Gradient
    # Obtain species KDEs for simulating from an empirical DCA-1 gradient
    
    if(length(gradientOrigDCA) != nrow(origAbundData)){
        stop("DCA and abundance data seem to not agree on number of samples/sites")
        }
    
    if(reportTests){
        # Sampling is not even across the DCA 1 gradient, as shown here:
        gradientHist <- hist(
            gradientOrigDCA, 
            breaks = nBreaksGradientHist, 
            plot = reportTests
            ,main = "Number of Samples in each Binned Segment of Gradient"
            ,xlab = "Gradient (DCA-1 Values)"
            )
    }else{
        gradientHist <- hist(
            gradientOrigDCA, 
            breaks = nBreaksGradientHist, 
            plot = reportTests
            #,main = "Number of Samples in each Binned Segment of Gradient"
            #,xlab = "Gradient (DCA-1 Values)"
            )
        }
    
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
    
    if(reportTests){
        # Some unneccessary plots 
        #heatmap(gradientOrigDCA, as.matrix(origAbundData[,1]))
        
        # how does abundance vary with DCA-1 for Epistominella pacifica?
        plot(gradientOrigDCA, 
            origAbundData$EpistominellaPacifica,
            xlab = "Original DCA Gradient", 
            ylab = "Single-Site Absolute Abundance of E. pacifica"
            )

        # how does abundance vary with DCA-1 for Lagena Spp.
        plot(gradientOrigDCA, 
            origAbundData$LagenaSpp,
            xlab = "Original DCA Gradient", 
            ylab = "Single-Site Absolute Abundance of Lagena Spp."
            )
        
                
        # how does RELATIVE abundance vary with DCA-1 for Epistominella pacifica?
        plot(gradientOrigDCA, 
            origAbundData_RelScale[,"EpistominellaPacifica"] 
                 / modeledSiteAbundance,
            xlab = "Original DCA Gradient", 
            ylab = "Single-Site Relative Abundance of E. pacifica"
            )
        
        # how does RELATIVE abundance vary with DCA-1 for Epistominella pacifica?
        plot(gradientOrigDCA, 
            origAbundData_RelScale[,"LagenaSpp"] 
                 / modeledSiteAbundance,
            xlab = "Original DCA Gradient", 
            ylab = "Single-Site Relative Abundance of Lagena Spp."
            )
        
        }
    
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
    
    # more plotting....
    if(reportTests){
        hist(gradientRepAbund$EpistominellaPacifica, 
            xlab = "Gradient (DCA 1 Score)", 
            breaks = histBreaks,
            main = paste0(
                "Bin-Summed Abundance (out of a ", modeledSiteAbundance,
                ") for E. pacifica\nSummed Across Samples in each Bin"
                )
            )
        
        hist(gradientRepAbund$LagenaSpp, 
            xlab = "Gradient (DCA 1 Score)", 
            breaks = histBreaks,
            main = paste0(
                "Expected Abundance (out of a ", modeledSiteAbundance,
                ") for Lagena Spp.\nSummed Across Samples in each Bin"
                )
            )
        }
    
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
    
    if(reportTests){
        hist(gradientRepSampAbund$EpistominellaPacifica, 
             xlab = "Gradient (DCA 1 Score)", breaks = histBreaks,
             main = paste0("Expected Per-", modeledSiteAbundance,
                "-Specimen-Sample Abundance for E. pacifica")
            )
        
        hist(gradientRepSampAbund$LagenaSpp, 
             xlab = "Gradient (DCA 1 Score)", breaks = histBreaks,
             main = paste0("Expected Per-", modeledSiteAbundance,
                "-Specimen-Sample Abundance for Lagena Spp.")
            )
        }
    
    # Now fit the KDE to each species' abundance curve
    kdePerTaxa <- lapply(gradientRepSampAbund, density, 
        bw = diff(gradientHist$mids)[1] )

    # Get maximum and min abundance for each species 
    # from the sampling normalized gradientSamplingCorrAbund 
    
    # Get mean relative abundance as well. 
    # Assign max/min and mean relative abundance to
        # elements stored within the list object from
        # the kernal density estimation function

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
        modeledSiteAbundance,
        reportTests = reportTests){
        
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
        modeledSiteAbundance = modeledSiteAbundance,
        reportTests = reportTests
        )
    
    if(reportTests){
        #Let's plot the KDEs for all taxa
        plot(kdeRescaled[[1]],
             main = "Scaled KDEs of Taxon Abundance",
             ylab = "Expected Relative Abundance",
             xlab = "Gradient (DCA1)",
             ylim = c(0, 0.6),
             #log = "y",
             xlim = c(min(gradientOrigDCA), max(gradientOrigDCA))
            )
        for(i in 2:length(kdeRescaled)){
            lines(kdeRescaled[[i]], col = i)
            }
        
        #Examine KDE for just *Epistominella pacifica*
        whichEP <- which(names(kdeRescaled) == "EpistominellaPacifica")
        plot(kdeRescaled[[whichEP]],
             main = "Scaled KDEs for Epistominella pacifica",
             ylab = "Expected Relative Abundance",
             xlab = "Gradient (DCA1)",
             lwd = 2,
             #log = "y",
             xlim = c(min(gradientOrigDCA), max(gradientOrigDCA))
            )
        
        whichLagSpp <- which(names(kdeRescaled) == "LagenaSpp")
        plot(kdeRescaled[[whichLagSpp]],
             main = "Scaled KDEs for Lagena Spp",
             ylab = "Expected Relative Abundance",
             xlab = "Gradient (DCA1)",
             lwd = 2,
             #log = "y",
             xlim = c(min(gradientOrigDCA), max(gradientOrigDCA))
            )
        }
    
    return(kdeRescaled)
    }

