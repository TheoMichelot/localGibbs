
#' Generate track from local Gibbs model
#' 
#' @param nbObs Number of locations to simulate
#' @param beta Parameters of the RSF
#' @param allr Vector of availability radii for simulation
#' @param covlist List of covariate rasters
#' @param xy0 Initial location (defaults to middle of map)
#' @param norm Logical. If TRUE, a normal transition density is used (with 
#' variance allr[1]^2).
#' 
#' @return Matrix of locations
#' 
#' @export
simLG <- function(nbObs, beta, allr, covlist, xy0=NULL, norm=FALSE)
{
    if(length(beta)!=length(covlist))
        stop("'beta' and 'covlist' should be the same length")
    if(length(allr)!=1 & length(allr)!=nbObs)
        stop("'allr' should be of length nbObs or 1")
    
    # store covariates in array
    ncov <- length(covlist)
    covvec <- NULL
    for(i in 1:ncov) {
        covmat <- raster::as.matrix(covlist[[i]])
        covvec <- c(covvec,covmat)
    }
    cov <- array(covvec, dim=c(dim(covlist[[1]])[1:2],ncov))
    
    # raster limits and resolution
    lim <- as.vector(extent(covlist[[1]]))
    res <- c(xres(covlist[[1]]), yres(covlist[[1]]))
    
    # xy0 default to middle point
    if(is.null(xy0))
        xy0 <- c((lim[1]+lim[2])/2, (lim[3]+lim[4])/2)
    
    if(length(allr)==1)
        allr <- rep(allr, nbObs)
    
    xy <- simLG_rcpp(nbObs=nbObs, beta=beta, allr=allr, cov=cov, xy0=xy0, 
                     lim=lim, res=res, norm=ifelse(norm,1,0))
    return(xy)
}
