
#' Categorical to dummy covariates
#' 
#' @param rast Raster of categorical covariate
#' 
#' @return List of dummy covariate rasters
#' 
#' @export
catToDummy <- function(rast)
{
    mat <- raster::as.matrix(rast)
    vals <- sort(unique(c(mat)))
    ncov <- length(vals)
    
    if(ncov>30)
        stop("Are you sure this is a categorical covariate? There are more than 30 different values.")
    
    # store covariates in array and list
    covlist <- list()
    for(i in 1:ncov) {
        h <- vals[i]
        covmat <- matrix(0, nrow(mat), ncol(mat))
        ind <- which(mat==h, arr.ind=TRUE)
        covmat[ind] <- 1
        
        covlist[[i]] <- rast # to have the right format
        values(covlist[[i]]) <- c(t(covmat)) # to have the right values
    }
    
    return(covlist)
}
