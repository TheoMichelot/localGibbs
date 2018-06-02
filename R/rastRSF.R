
#' Calculate raster of RSF
#' 
#' @param beta Coefficients of the RSF
#' @param covlist List of covariate rasters
#' 
#' @return Raster of the RSF
#' @export
rastRSF <- function(beta, covlist)
{
    if(length(beta)!=length(covlist))
        stop("'beta' and 'covlist' sould have the same length")
    
    rast <- 0
    for(i in 1:length(covlist))
        rast <- rast + beta[i]*covlist[[i]]
    
    return(exp(rast))
}
