#' Plot Kernel Density Estimates of Species Abundance Across a Focal Gradient
#' 
#' The

#' @details
#' The

# @inheritParams

#' @param speciesKDEs The

#' @param fullGradientRange The

#' @param xlim The

#' @param ylim The

#' @param logY The

#' @return
#' The

# @aliases

# @seealso

# @references

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
             ylab = "Scaled Expected Abundance If Present",
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
             xlab = "Gradient (DCA-1)",
             xlim = xlim,
             ylim = ylim,
             lwd = 2
             #,log = "y"
             )
        }

    for(i in 2:length(speciesKDEs)){
        lines(speciesKDEs[[i]], 
            col = i, lwd = 2
            )
        }

    }