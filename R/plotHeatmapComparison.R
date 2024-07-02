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
            key.axes = NULL,
        
            palette = "plasma",
    
            xlab = "x", 
            ylab = "y",
            main = "main title",
        
            gradientKeyLabel = "color gradient key",
            mtext_line = 3,
            margins = c(5,6,5,5)
            ){

    ########################
    par(mar = margins)

    if(is.null(contourLevels)){
        contourLevels <- pretty(zlim, nlevels)
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
        
    filled.contour(
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
            hcl.colors(n, palette = palette), 
        
        plot.axes = { 
            
            if(!is.null(contourLineLevels)){
                contour(
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
            
            axis(1, at=xtickAt, label=xtick); 
            axis(2, at=ytickAt, label=ytick);
            
            if(!is.null(additionalFeatures)){
                additionalFeatures
                }
            }

        # this code seems to stop contours from being filled...
        #,if(!is.null(key.axes)){
        #    key.axes = key.axes
        #    }

                )

    mtext(
        side = 4, cex = 1, 
        text = gradientKeyLabel,
        line = mtext_line
        )
    }

