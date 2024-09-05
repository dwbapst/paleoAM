
# Loosely based on vegan:::predict.decorana                
    
#  an example of stuff I don't need from 
            # vegan:::predict.decorana
            
# u <- sweep(x = oldDCA$rproj[ , 1, drop = FALSE], 
#           MARGIN = 2, 
#           STATS = oldDCA$origin[1], 
#           FUN = "-")
#
# v <- sweep(x = oldDCA$cproj[ , 1, drop = FALSE], 
#           MARGIN = 2, 
#           STATS = oldDCA$origin[1], 
#           FUN = "-") 
# rs <- oldDCA$aidot/sum(oldDCA$aidot)                
# u <- oldDCA$rproj[ , 1] - oldDCA$origin[1]

predictDCA_sites_dca1 <- function (
            oldDCA, 
            newSample
            ){
    cs <- oldDCA$adotj/sum(oldDCA$aidot)
    v <- oldDCA$cproj[ , 1] - oldDCA$origin[1]
    v <- v *  sqrt(cs)    
    #
    proj <- as.matrix(newSample)
    if (!is.null(oldDCA$v)){
        proj <- sweep(proj, 2, oldDCA$v, "*")
        }
    rs <- rowSums(proj)
    proj <- (proj - outer(rs, cs))/sqrt(outer(rs, cs))
    out <- sweep(proj %*% v, 1, sqrt(rs), "/")
    return(out)
    }