#' Plots a Heatmap Comparison
#' 
#' This function is a complex wrapper of the functions \code{contour} and \code{filled.contour} which allows a user to combine colored contours for the surface of a third variables plotted across a two-dimensional space, with contour lines of the third variable's surface also plotted in that same two-dimensional space. This is often a common plotting need for this package, and thus is included here.

#' @details
#' The function \code{filled.contour} doesn't easy allow for sequential modifications, like adding additional contour lines to an existing contour plot, and so this function simplifies having to write the second \code{contour} plot as an argument for \code{plot.axes} in \code{filled.contour}.

# @inheritParams

#' @param x,y The horizontal (\code{x}) and vertical (\code{y}) variables that determine the two-dimensional space within which the third variable (\code{z}) is plotted as a surface.

#' @param z The values of the third variable (\code{z}) that will be used to define the plotted surface for contours, given as a matrix with the same number of rows as the length of \code{x}, and the same number of columns as the length of \code{y}.

#' @param xlim,ylim,zlim These are two-element vectors giving the minimum and maximum limits for the horizontal variable (\code{x}), vertical variable (\code{y}), and the variable defining the surface plotted (\code{z}) within the two dimensional space defined by \code{x} and \code{y}.

#' @param xlog,ylog Should the \code{x} or \code{y} axes be displayed with log-scaling?

#' @param xtick,ytick Vectors that give the positions of the tick-marks for \code{x} and \code{y} axes.

#' @param contourLevels A vector of values at which to put the breaks between the color-filled contours for \code{z}, Also determines the
#' different levels show on the color-gradient key shown to the side of the contour plot.

#' @param nlevels The number of different color levels to use, if \code{contourLevels} is not defined.

#' @param contourLineLevels A vector of values at which to put the distinct contour lines for \code{z}. This must be defined for contour lines to be plotted.

#' @param contour.lwd The thickness of plotted contour lines.

#' @param additionalFeatures Additional features to add to the contour space, such as 

# @param key.axes ???

#' @param palette The color palette to use for \code{filled.contour}. By default, the palette \code{"plasma"} is used.

#' @param xlab,ylab The labels for the \code{x} and \code{y} axes.

#' @param main The plot's main title.

#' @param gradientKeyLabel The optional label text for the color gradient key for \code{z}, shown to the right of the main contour plot. This label will be shown to the right of the key.

#' @param mtext_line The distant in the margin away from the key at which the \code{gradientKeyLabel} is displayed. The default value is 3.

#' @param margins The size of the margins for the result plot. The default configuration gives some extra room on the left-hand size.


#' @return
#' This function returns nothing at all as output. It just makes a plot.

# @aliases

# @seealso

# @references

# @examples

#' @name plotHeatmapComparison
#' @rdname plotHeatmapComparison
#' @export
plotHeatmapComparison <- function(
    
            x, y, z,  
    
            xlim = range(x, finite = TRUE), 
            ylim = range(y, finite = TRUE), 
            zlim = range(z, finite = TRUE),
    
            xlog = FALSE, 
            ylog = FALSE,
    
            xtick = pretty(x),
                #c(0.02, 0.05, 0.1, 0.2, 0.3, 
                #0.4, 0.5, 0.7, 1, 2, 5, 10),
            ytick = pretty(y),
                #c(0.02, 0.05, 0.1, 0.2, 0.3, 
                #0.4, 0.5, 0.7, 0.9),

            contourLevels = NULL,
            nlevels = 10,
    
            contourLineLevels = NULL,
            contour.lwd = 2,
    
            additionalFeatures = NULL,
            #key.axes = NULL,
        
            palette = "plasma",
    
            xlab = "x", 
            ylab = "y",
            main = "main title",
        
            gradientKeyLabel = "color gradient key",
            mtext_line = 3,
            margins = c(5,6,5,5)
            ){

    ########################
    # setting up par to reset on exit 
        # as instructed by CRAN
    oldpar <- par(no.readonly = TRUE) # code line i
    on.exit(par(oldpar)) # code line i + 1
    
    ########################    
    graphics::par(mar = margins)

    if(is.null(contourLevels)){
        contourLevels <- base::pretty(zlim, nlevels)
        }
    
    xtickAt <- xtick
    ytickAt <- ytick
    
    if(xlog){
        x <- log(x, 10)
        xtickAt <- log(xtick, 10)
        xlim <- log(xlim, 10)
        }
    if(ylog){
        y <- log(y, 10)
        ytickAt <- log(ytick, 10)
        ylim <- log(ylim, 10)
        }
        
    graphics::filled.contour(
        x = x,
        y = y,
        z = z,
        
        xlim = xlim,
        ylim = ylim,
        zlim = zlim,
        
        xlab = xlab,
        ylab = ylab,
        main = main,
        
        # control the breaking up of contours for the z surface
        levels = contourLevels,
        nlevels = nlevels,
        
        color.palette = function(n) 
            grDevices::hcl.colors(n, palette = palette), 
        
        plot.axes = { 
            
            if(!is.null(contourLineLevels)){
                graphics::contour(
                    x = x,
                    y = y,
                    z = z,
                    
                    xlim = xlim,
                    ylim = ylim,
                    zlim = zlim,
                    
                    levels = contourLineLevels, 
                    lwd = contour.lwd,
                    
                    drawlabels = FALSE, 
                    axes = FALSE, 
                    frame.plot = FALSE, 
                    add = TRUE
                    )
                };
            
            graphics::axis(1, at=xtickAt, label=xtick); 
            graphics::axis(2, at=ytickAt, label=ytick);
            
            if(!is.null(additionalFeatures)){
                additionalFeatures
                }
            }

        # this code seems to stop contours from being filled...
        #,if(!is.null(key.axes)){
        #    key.axes = key.axes
        #    }

                )

    graphics::mtext(
        side = 4, cex = 1, 
        text = gradientKeyLabel,
        line = mtext_line
        )
    }

