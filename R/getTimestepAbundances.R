#' Simulate Fossil Assemblages with Abundances at each Time-Step
#'
#' Given a set of KDEs fit to species abundance and models of species occurrance relative to an environmental gradient, and given a sequence of gradient values, and a number of specimens to sample at each time-step, obtains a matrix containing abundances for species as a series of simulated assemblages. 

#' @details
#' \code{getTimestepAbundances} represents simulating the original biotic community that was present at some given point in time, which is not the same thing as a fossil assemblage that might be collected from sediments today as finite samples. That is covered by feeding the output from this function to \code{sampleFossilAssemblageSeries}.

# @inheritParams simulateFossilAssemblageSeries

#' @param kdeRescaled The list of modeled KDEs for species abundance, output from \code{\link{getSpeciesSpecificRescaledKDE}}.

#' @param probSpeciesOccur The output from \code{\link{getProbOccViaPresAbs}}

#' @param gradientValues A vector of gradient values to simulate over. A separate 'true' assemblage / community will be simulated for each value in the respective vector.

#' @param specimensPerTimestep The number of specimens returned in a given time-step by \code{getTimestepAbundances}, usually set to an unrealistically high number to represent the true 'unsampled' fossil assemblage.

#' @return
#' A matrix containing abundances for species as a series of simulated assemblages.

# @seealso

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
        approx(x = kdeEst$x, y = kdeEst$y, xout = unqGradient)$y               
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
    # 05-19-21
    # only certain species can be sampled at some point along the gradient
        # use proportion of samples in each 0.2 bin as a rough approximation    
    #OLD: speciesPresent_List <- lapply(unqGradient, probSpeciesOccur)
    speciesPresent_List <- probSpeciesOccur(unqGradient) 
    #that produce a list that is nspecies long, 
        # with each species-specific element 
        # containing values for each unique gradient value
    
    # test
    if(length(speciesPresent_List) != nSpecies){
        stop("not getting right number of species from probSpeciesOccur")
        }
    
    # need to make it a list with elements for each unique gradient value...
    speciesPresent_List <- lapply(1:nUnqGradient, function(x) 
        lapply(speciesPresent_List, function(y) y[[x]]))
    # and now make it a list for each timestep, using matchUnqGradient
    speciesPresent_List <- speciesPresent_List[matchUnqGradient]
    # now it is a list that has elements for each timestep, 
        # with each element being nSpecies long
    # simplify
    speciesPresent_List <- lapply(speciesPresent_List, unlist)
    
    # test that its the right length
    if(length(speciesPresent_List) != nTimeSteps){
        stop("not getting right number of values for speciesPresent_List")
        }
    
    # now stochastically determine if species are present or not
    # first generate a large matrix of numbers pulled from a uniform distribution
    uniformDistNumbers <- runif(n = nSpecies * nTimeSteps, 
         min = 0, max = 1)
    uniformDistNumbers <- matrix(uniformDistNumbers, 
         nTimeSteps, nSpecies)
    
    # sample from a uniform distribution (0 -> 1)
        # as a way of getting stochastic presence/absence
    speciesPresent_List <- lapply(1:nTimeSteps, function(i){
            probs <- speciesPresent_List[[i]]
            uniformDraw <- uniformDistNumbers[i,]
            uniformDraw <= probs
            })

    # now simulate actual specimen abundances for each time step    
    timestepAbundances <- simulateTimestepAbundances(
        specimensPerTimestep = specimensPerTimestep, 
        nSpecies = nSpecies, 
        nTimeSteps = nTimeSteps, 
        speciesPresent_List = speciesPresent_List, 
        expAbundFromKDE_List = expAbundFromKDE_List
        )
    
    return(timestepAbundances)
    }



simulateTimestepAbundances <- function(
            specimensPerTimestep, 
            nSpecies, 
            nTimeSteps, 
            speciesPresent_List, 
            expAbundFromKDE_List){
    
    # make empty abundance matrix
    timestepAbundances <- matrix( 0, 
          nrow = nTimeSteps, ncol = nSpecies)
    
    for(i in 1:nTimeSteps){
        # figure out relative expected frequency 
        # of each species at each point in time
        expAbundFromKDE <- expAbundFromKDE_List[[i]]    
        # conditional on IF they were sampled    
        speciesPresent <- speciesPresent_List[[i]]        
        # retain only present species
        # set expected abundance of all other species to 0
        expAbundFromKDE[!speciesPresent] <- 0
        # turn into expected relative abundances
        expRelativeAbundances <- expAbundFromKDE/sum(expAbundFromKDE)
        
        # sample "specimensPerTimestep" fossil specimens for each timestep
        # We treat each timestep as having a fixed non-stochastic number of individuals 
        # sampled from that community (specimensPerTimestep) 
        
        species <- (1:nSpecies)[expRelativeAbundances > 0]
        species <- as.integer(species)
        expRelativeAbundances <- expRelativeAbundances[expRelativeAbundances > 0]
            
        fossilSamples <- sample(size = specimensPerTimestep, 
                x = species, replace = TRUE, 
                prob = expRelativeAbundances)
        
        # count how many of each species were buried
        fossilCounts <- tabulate(fossilSamples, nbins = nSpecies)
        timestepAbundances[i,] <- fossilCounts
        }

    return(timestepAbundances)
    }
