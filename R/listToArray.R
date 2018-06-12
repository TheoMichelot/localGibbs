
#' Transform list of rasters to array
#' 
#' @param covlist List of covariate rasters
#' 
#' @return Array, for which each layer corresponds to a covariate
#' 
#' @export
listToArray <- function(covlist)
{
    ncov <- length(covlist)
    covvec <- NULL
    for(i in 1:ncov) {
        covmat <- raster::as.matrix(covlist[[i]])
        covvec <- c(covvec,covmat)
    }
    covarray <- array(covvec, dim=c(dim(covlist[[1]])[1:2],ncov))
    
    return(covarray)
}
