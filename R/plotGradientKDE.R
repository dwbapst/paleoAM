#' Plot Kernel Density Estimates of Species Abundance Across a Focal Gradient
#' 
#' This function plots each rescaled KDE fit to each specific-specific 
#' rise-and-fall in abundance across some ecological gradient variable.

#' @details
#' In many ways, this is an attempt to create empirical versions of the the
#' hypothetical figures portraying abundance response curves 
#' relative to an environmental gradient in Patzkowsky & Holland (2012).
#' 
#' The ecological gradient variable is often an environmental characteristic, 
#' such as depth, oxygenation, altitude, precipitation, 
#' but this is not necessarily so. 

# @inheritParams

#' @param speciesKDEs A list of rescaled-KDE data, where each element is a 
#' different species, such as that output by \code{\link{getSpeciesSpecificRescaledKDE}}.

#' @param fullGradientRange The minimum and maximum value of the ecological 
#' gradient variable at which ecological assemblage data was observed. 
#' If \code{xlim} isn't given, this defines the horizontal axis limits for resulting plot.

#' @param xlim,ylim Vectors of two elements, defining the minimum and maximum 
#' values for the horizontal (x) axis and vertical (y) axis, respectively. 
#' The default for \code{xlim} is \code{NULL} and only needs to be defined 
#' if different axis limits than \code{fullGradientRange} is desired. 
#' The default for \code{ylim} is \code{c(0,1)} which likely leaves considerable
#'  empty white space above the KDEs, which can be reduced by adjusting this argument.

#' @param logY Should the vertical axis (the relative height of rescaled KDEs) 
#' be portrayed with logarithmic scaling?

#' @return
#' Nothing is returned, just a plot is made.

# @aliases

#' @seealso
#' This function mainly exists to look at the output from 
#' \code{\link{getSpeciesSpecificRescaledKDE}} for a fossil assemblage.

#' @references
#' Patzkowsky, M.E. and Holland, S.M., 2012. \emph{Stratigraphic Paleobiology: 
#' Understanding the Distribution of Fossil Taxa in Time and Space.}
#' University of Chicago Press. 259 pages.

# @examples

#' @name plotGradientKDE
#' @rdname plotGradientKDE
#' @export
plotGradientKDE <- function(
        speciesKDEs,
        fullGradientRange,
        xlim = NULL,
        ylim = c(0, 1),
        logY = FALSE
        ){

    if(is.null(xlim)){
        xlim <- fullGradientRange
        }
    
    if(logY){
        plot(speciesKDEs[[1]],
             #main = "Scaled KDEs of Abundance for Every Species in the Dataset",
             main = "",
             ylab = "Log Scaled Expected Abundance\n(If Present)",
             xlab = "Gradient (DCA-1)",
             xlim = xlim,
             ylim = ylim,
             lwd = 2
             ,log = "y"
             )
    }else{
        plot(speciesKDEs[[1]],
             #main = "Scaled KDEs of Abundance for Every Species in the Dataset",
             main = "",
             ylab = "Scaled Expected Abundance If Present",
             # ylab = "Scaled Expected Abundance\n(If Present)",
             xlab = "Gradient (DCA-1)",
             xlim = xlim,
             ylim = ylim,
             lwd = 2
             #,log = "y"
             )
        }

    for(i in 2:length(speciesKDEs)){
        graphics::lines(speciesKDEs[[i]], 
            col = i, lwd = 2
            )
        }

}
