
#' Plot raster with ggplot
#' 
#' @param rast Raster object
#' @param xy Optional matrix of locations
#' @param norm Logical. Should the raster be normalized to integrate to 1?
#' @param log Logical. Should the raster be plotted on the log scale?
#' @param name Name of the scale
#' 
#' @importFrom ggplot2 ggplot geom_raster coord_equal geom_point geom_path aes_string
#' @importFrom viridis scale_fill_viridis
#' @importFrom raster xres yres values 
#' @importFrom sp coordinates
#' 
#' @export
plotRaster <- function(rast, xy=NULL, norm=FALSE, log=FALSE, name="")
{
    covmap <- data.frame(coordinates(rast),val=values(rast))
    
    if(norm) {
        s <- sum(covmap$val) * xres(rast) * yres(rast)
        covmap$val <- covmap$val/s
    }
    if(log) {
        covmap$val <- log(covmap$val)
    }
    
    p <- ggplot(covmap, aes_string(x="x",y="y")) + geom_raster(aes_string(fill="val")) +
        coord_equal() + scale_fill_viridis(name=name)
    
    if(!is.null(xy)) {
        xydf <- data.frame(x=xy[,1], y=xy[,2])
        p <- p + geom_point(aes_string(x="x",y="y"), data=xydf, size=0.3) +
            geom_path(aes_string(x="x",y="y"), data=xydf, size=0.5)
    }
    
    return(p)
}
