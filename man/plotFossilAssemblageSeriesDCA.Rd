% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotFossilAssemblageSeriesDCA.R
\name{plotFossilAssemblageSeriesDCA}
\alias{plotFossilAssemblageSeriesDCA}
\title{Plot the Recovered DCA Values of a Fossil Assemblage Series}
\usage{
plotFossilAssemblageSeriesDCA(
  ...,
  colSimGenerating = "black",
  colSimRecovered = "navy"
)
}
\arguments{
\item{...}{This function takes either the output from 
\code{\link{simulateFossilAssemblageSeries}} or requires specifying 
the three arguments \code{simTimeVar}, a data-frame of time-steps, 
sedimentary width and gradient values for a time-series simulation); 
\code{gradientRecovered}, the recovered gradient values; and
\code{sampleAge}, the age of individual samples.}

\item{colSimGenerating}{What color should be used for the generating 
("true") gradient values?}

\item{colSimRecovered}{What color should be used for the recovered 
gradient values calculated from the simulated data?}
}
\value{
Returns nothing at all. Just a plot. That's all!
}
\description{
Makes a plot of the simulated generating gradient over time along with the
recovered gradient values in the same plot for a given simulated time series 
of fossil assemblages, particularly those from
\code{\link{simulateFossilAssemblageSeries}}.
}
\details{
The function will generally only be run on the output 
from \code{\link{simulateFossilAssemblageSeries}}, although the function is 
written so that the necessary elements can be provided separately.
}
