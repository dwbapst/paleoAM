test_that("quick example of simulateFossilAssemblageSeries with bioturbation works", {

data(gulfOfAlaska)

alaskaKDEs <- getSpeciesSpecificRescaledKDE(
    gradientOrigDCA = DCA1_GOA, 
    origAbundData = abundData_GOA, 
    abundanceFloorRatio = 0.5, 
    nBreaksGradientHist = 20, 
    modeledSiteAbundance = 10000
    )
    
alaskaProbOccur <- getProbOccViaPresAbs(
   gradientOrigDCA = DCA1_GOA, 
   origAbundData = abundData_GOA
   )

# Run the simulation of fossil assemblages
    # simulateFossilAssemblageSeries has lots of arguments...
    # below they are broken up into groups, seperate by #
    # matches scenarios from fig 13 of Belanger & Bapst
    
fossilSeriesOut <- simulateFossilAssemblageSeries(
      # inputs
      kdeRescaled = alaskaKDEs,
      probSpeciesOccur = alaskaProbOccur,
      origAbundData = abundData_GOA,
      fullGradientRange = c(min(DCA1_GOA), max(DCA1_GOA)),
      
      # let's make it relatively mild event 
        # with a long transition 
      eventChangeScale = 0.5,
      bgGradientValue = -1,
      transitionDurationRatio = 1,
       
      # don't need to define eventSampleWidthRatio 
        # - only need to define three of eventSampleWidthRatio, 
        # sampleWidth, eventDuration, sedRatePerTimestep
      sampleWidth = 3,
      eventDuration = 100, 
      sedRatePerTimestep = 0.1,
      
      # sample every third sample-width worth of core
      samplingCompleteness = 1/3,
      # no bioturbation 
      bioturbDepthRatio = 10/3,
      bioturbIntensity = 1,
           
      nEvents = 1,
      nSpecimens = 100,
      # let's plot it
      plot = TRUE
      )

  expect_type(fossilSeriesOut,type = "list") 

})
