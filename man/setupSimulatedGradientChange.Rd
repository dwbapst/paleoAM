% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setupSimulatedGradientChange.R
\name{setupSimulatedGradientChange}
\alias{setupSimulatedGradientChange}
\title{Create a Stochastic Time-Series of Gradient Change For Use in Simulating Assemblage Change}
\usage{
setupSimulatedGradientChange(
  nEvents,
  peakGradientValue,
  bgGradientValue,
  bgDurationRange,
  transitionDuration,
  eventDuration,
  halfGradientOnly = FALSE,
  includeInitialBackgroundPhase = TRUE,
  plot = FALSE
)
}
\arguments{
\item{nEvents}{Number of events to occur in a simulated sequence 
of gradient change.}

\item{peakGradientValue}{The gradient value at the 'peak' for 
an event that represents an excursion on that environmental gradient.}

\item{bgGradientValue}{The gradient value expected during 
background intervals during which no notable excursion is 
occurring on that environmental gradient.}

\item{bgDurationRange}{A vector of two values, representing the 
minimum and maximum duration (in time units) of a 
background interval between successive events.}

\item{transitionDuration}{How long the transition between peak and 
background intervals should be. The longer this interval, 
the more chances of an assemblage being sampled that 
represents transitional gradient values.}

\item{eventDuration}{The duration (in time-units) of a 
simulated event during which the environmental gradient 
is at an excursion 'peak' level.}

\item{halfGradientOnly}{Whether to simulate only half of 
a background-event sequence, either beginning or terminating 
the simulation at the peak value. 
Only a single event can be simulated, so \code{nEvents} must be 1. 
The default is \code{FALSE} which signals to not simulated a half-gradient.}

\item{includeInitialBackgroundPhase}{A logical indicating whether 
to include a lengthy background phase, for use in calibrating a simulation. 
This function is mainly for diagnostic purposes 
and may be removed in future updates.}

\item{plot}{Should the simulated gradient be shown as a plot?}
}
\value{
A list with five components: 
\code{simGradient}, a data frame giving the change in gradient values over time; 
\code{approxGradientSeriesFunction}, the simulated gradient given as an interpolated function;
\code{eventStartEndTimes}, a vector of when each event and its preceding transition begin in time-units;
\code{eventPhaseStartTimes}, a vector of when each new event phase (at the peak gradient value) begin in time-units;
and \code{backgroundStartEnd}, a value indicating the time-step when the beginning background interval ends.
}
\description{
Given a series of inputs, simulates a sequence of gradient change 
against time for use in testing how environmental change alters 
the recovered sequence of fossil assemblages.
}
\details{
This function is rather complicated and was written at a time 
when it was envisioned that simulations would involve time series 
of many repeated events with varying background intervals 
between them, rather than simulated sequences having only one event. 
In practice, use of paleoAM has tended to find the latter to be more useful.
}
