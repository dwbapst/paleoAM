#' Simulate Fossil Assemblages with Abundances at each Time-Step
#'
#' Given a set of KDEs fit to species abundance and models of species occurrence relative to an environmental gradient, and given a sequence of gradient values, and a number of specimens to sample at each time-step, obtains a matrix containing abundances for species as a series of simulated assemblages. 

#' @details
#' \code{getTimestepAbundances} represents simulating the original biotic community that was present at some given point in time, which is not the same thing as a fossil assemblage that might be collected from sediments today as finite samples. That is covered by feeding the output from this function to \code{sampleFossilSeries}.
#' 
#' Thus, this function is generally run before running \code{\link{sampleFossilSeries}}, 
#' however most users will likely never run either function, 
#' instead running \code{\link{simulateFossilAssemblageSeries}}. 

# @inheritParams simulateFossilAssemblageSeries

#' @param kdeRescaled The list of modeled KDEs for species abundance, output from \code{\link{getSpeciesSpecificRescaledKDE}}.

#' @param probSpeciesOccur The output from \code{\link{getProbOccViaPresAbs}}

#' @param gradientValues A vector of gradient values to simulate over. A separate 'true' assemblage / community will be simulated for each value in the respective vector.

#' @param specimensPerTimestep The number of specimens returned in a given time-step by \code{getTimestepAbundances}, usually set to an unrealistically high number to represent the true 'unsampled' fossil assemblage.

#' @return
#' A matrix containing abundances for species as a series of simulated assemblages.

#' @seealso
#' This function is generally run before running \code{\link{sampleFossilSeries}}. 
#' Most users will likely never run either function, instead running \code{\link{simulateFossilAssemblageSeries}}.

# @references

# @examples

#' @name getTimestepAbundances
#' @rdname getTimestepAbundances
#' @export
getTimestepAbundances <- function(
              kdeRescaled, 
              probSpeciesOccur, 
              gradientValues, 
              specimensPerTimestep
              ){
    
    nSpecies <- length(kdeRescaled)
    nTimeSteps <- length(gradientValues)
    # how many unique gradient values are there?
    unqGradient <- sort(unique(gradientValues))
    nUnqGradient <- length(unqGradient)
    # match unique gradient values to the simulated gradient curve
    #matchUnqGradient <- lapply(gradientValues, function(x) 
    #    which(unqGradient == x)[1])  
    matchUnqGradient <- match(gradientValues, unqGradient)
    matchUnqGradient <- unlist(matchUnqGradient)
    
    # Generate KDEs before running simulation
    # Use approx function for estimating KDEs from gradient values
    # function for getting exp abundance from rescaled KDEs
    expAbundFromKDE_List <- lapply(kdeRescaled, function(kdeEst)
        stats::approx(x = kdeEst$x, y = kdeEst$y, xout = unqGradient)$y               
        )
    
    # test
    if(length(expAbundFromKDE_List) != nSpecies){
        stop("not getting right number of species from kdeRescaled")
        }
    # what is the structure of 
    if(length(expAbundFromKDE_List[[1]]) < 1){
        stop("expAbundFromKDE_List elements do not *any* contain values")
        }
    
    # this gives us a list that has elements for each species,
        # each composed of values for unique gradient values
    # flip it: now a list of unique gradient values with species values    
    expAbundFromKDE_List <- lapply(1:nUnqGradient, function(x)
        lapply(expAbundFromKDE_List, "[[", x))
    
    # using gradient matches, get expected abundances for each time step
    expAbundFromKDE_List <- expAbundFromKDE_List[matchUnqGradient]
    # simplify
    expAbundFromKDE_List <- lapply(expAbundFromKDE_List, unlist)        
    
    # get species occurrences for each time step based on occurrence data
    # 05-19-21: this is now unnecessary
    # only certain species can be sampled at some point along the gradient
        # use proportion of samples in each 0.2 bin as a rough approximation    
    #OLD: speciesPresent_Matrix <- lapply(unqGradient, probSpeciesOccur)
    speciesPresent_Matrix <- probSpeciesOccur(unqGradient) 
    #that produce a matrix that has nrow = nspecies, 
        # values across the gradient are each column
        # (ncol = length(unqGradient) (which is nUnqGradient)
    
    # test
    if(nrow(speciesPresent_Matrix) != nSpecies){
        stop("not getting right number of species from probSpeciesOccur")
        }
    if(ncol(speciesPresent_Matrix) != nUnqGradient){
        stop("not getting right number of gradient values from probSpeciesOccur")
        }
    
    # repeat each column for each unique gradient value
    speciesPresent_Matrix <- speciesPresent_Matrix[,matchUnqGradient]
    
    # 07-16-25: this is old
    # need to make it a list with elements for each unique gradient value...
        # speciesPresent_List <- lapply(1:nUnqGradient, function(x) 
            #     lapply(speciesPresent_List, function(y) y[[x]]))
        # and now make it a list for each timestep, using matchUnqGradient
        #speciesPresent_Matrix <- speciesPresent_List[matchUnqGradient]
        # now it is a list that has elements for each timestep, 
            # with each element being nSpecies long
        # simplify
        #speciesPresent_Matrix <- lapply(speciesPresent_List, unlist)
        
    # test that its the right length
    if(ncol(speciesPresent_Matrix) != nTimeSteps){
        stop("not getting right number of values for speciesPresent_Matrix")
        }
    
    # now stochastically determine if species are present or not
    # first generate a large matrix of numbers pulled from a uniform distribution
    uniformDistNumbers <- matrix(
        stats::runif(
            n = nSpecies * nTimeSteps, 
            min = 0, max = 1
            ), 
        nSpecies, nTimeSteps) # rows = species, cols = timesteps

    speciesPresent_Matrix <- uniformDistNumbers <= speciesPresent_Matrix    
        
    ## sample from a uniform distribution (0 -> 1)
    #    # as a way of getting stochastic presence/absence
    #speciesPresent_List <- lapply(1:nTimeSteps, 
    #    function(i){
    #        probs <- speciesPresent_List[[i]]
    #        uniformDraw <- uniformDistNumbers[i,]
    #        uniformDraw <= probs
    #        }
    #    )

    # do I need to convert speciesPresent_Matrix to speciesPresent_List?
        # maybe not...

    # now simulate actual specimen abundances for each time step    
    timestepAbundances <- simulateTimestepAbundances(
        specimensPerTimestep = specimensPerTimestep, 
        nSpecies = nSpecies, 
        nTimeSteps = nTimeSteps, 
        speciesPresent_Matrix = speciesPresent_Matrix, 
        expAbundFromKDE_List = expAbundFromKDE_List
        )
    
    return(timestepAbundances)
    }


# NOT EXPORTED
simulateTimestepAbundances <- function(
            specimensPerTimestep, 
            nSpecies, 
            nTimeSteps, 
            speciesPresent_Matrix, 
            expAbundFromKDE_List){
    
    # make empty abundance matrix
    timestepAbundances <- matrix( 0, 
          nrow = nTimeSteps, 
          ncol = nSpecies)
    
    for(i in 1:nTimeSteps){
        # figure out relative expected frequency 
        # of each species at each point in time
        expAbundFromKDE <- expAbundFromKDE_List[[i]]    
        # conditional on IF they were sampled    
        # speciesPresent <- speciesPresent_Matrix[,i]        
        # retain only present species
        # set expected abundance of all other species to 0
        expAbundFromKDE[!speciesPresent_Matrix[,i] ] <- 0
        # turn into expected relative abundances
        expRelativeAbundances <- expAbundFromKDE/sum(expAbundFromKDE)

        # sample "specimensPerTimestep" fossil specimens for each timestep
        # We treat each timestep as having a 
        # fixed non-stochastic number of individuals 
        # sampled from that community (specimensPerTimestep) 
        
        # 07-15-25
        # following code not needed, can replace sample+tabulate
        # use a call to rmultinom instead
        #
            #species <- which(expRelativeAbundances > 0)
            #if(length(species) > 1){
                #species <- as.integer(species)
            #    fossilSamples <- sample.int(
            #        n = length(species),
            #        size = specimensPerTimestep, 
            #        # x = species, # x = as.character(species), 
            #        replace = TRUE, 
            #        prob = expRelativeAbundances[species]
            #        )
             #   fossilSamples <- species[fossilSamples]
                # convert fossilSamples to integer
                  # this eats up an enormous amount of computational cycles -- why??
                # fossilSamples <- as.integer(fossilSamples)
            # }else{
                # if there is only one species
              #  fossilSamples <- rep(species, specimensPerTimestep)
              #  }
            
        
        # count how many of each species were buried
        #fossilCounts <- tabulate(fossilSamples, nbins = nSpecies)
        
        timestepAbundances[i,] <- rmultinom(n = 1, 
                                            size = specimensPerTimestep, 
                                            prob = expRelativeAbundances
                                            )[,1]
        }

    return(timestepAbundances)
    }
