simulateGradientQuantile <- function(
            quantileProbs = c(0.95),
            nSamplesSim,
            gradientValue,
            origAbundData,
            kdeRescaled,
            probSpeciesOccur,
            specimensPerTimestep,
            nSpecimensPicked
            ){
    
    # previously named as generateBackgroundOnlySimulationQuantile
    
    # generate simulations based only on the background gradient value
        # use these to estimate quantiles which we use to denote 'significance'
        # for the purpose of testing whether an event sample was significant different
        # from the background values
    
    # Figure out which samples will be taken and at what core depths / outcrop heights 
    # First, make a matrix with age, depth, gradient values
    simTimeVar <- data.frame(
        timestep = 1:(nSamplesSim*3),
        # reverse core depth so oldest time is at bottom (biggest depths)
          # 07-08-21: subtract one to start from 0
        sedColumnDepth = ((1:nSamplesSim) -1) * 3,
        gradientValue = gradientValue
        )
        
    # Simulate actual specimen abundances for each time step
    # using KDEs and occurrence probabilities from empirical data
    timestepAbundances <- getTimestepAbundances(
        gradientValues = simTimeVar$gradientValue,
        kdeRescaled = kdeRescaled,
        probSpeciesOccur = probSpeciesOccur,
        specimensPerTimestep = specimensPerTimestep
        )
    colnames(timestepAbundances) <- names(kdeRescaled)
        
    # Running the Simulation
    fossilSeries <- sampleFossilAssemblageSeries(    
        bioturbIntensity = 0, 
        bioturbZoneDepth = 0,  
        distBetweenSamples = 0, 
        sampleWidth = 3, 
        simTimeVar = simTimeVar,
        timestepAbundances = timestepAbundances,
        nSpecimensPicked = nSpecimensPicked
        )
        
    # get sample DCA1 scores
    ecologyOutList <- quantifyCommunityEcology(
        origAbundData = origAbundData,
        fossilSeries = fossilSeries,
        singularDCA = TRUE,
        inclusiveDCA = FALSE,
        rawDCA = FALSE
        )

    DCA1 <- ecologyOutList$scoreDCA1_singular
    quantileOut <- quantile(
        x = DCA1, probs = quantileProbs, na.rm = TRUE
        )
    return(quantileOut)
    }
