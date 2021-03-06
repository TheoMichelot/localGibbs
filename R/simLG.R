
#' Generate track from local Gibbs model
#' 
#' @param nbObs Number of locations to simulate
#' @param beta Parameters of the RSF
#' @param allr Vector of availability radii (or standard deviations if norm=TRUE)
#' @param covlist List of covariate rasters
#' @param xy0 Initial location (defaults to middle of map)
#' @param norm Logical. If TRUE, a normal transition density is used (with 
#' variance given by allr^2).
#' @param npts Number of potential endpoints to sample at each time step
#' 
#' @return Matrix of locations
#' 
#' @export
simLG <- function(nbObs, beta, allr, covlist, xy0=NULL, norm=FALSE, npts=100)
{
    if(length(beta)!=length(covlist))
        stop("'beta' and 'covlist' should be the same length")
    
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
    
    if(length(allr)==1) {
        allr <- rep(allr, nbObs-1)
    } else if(length(allr)>=nbObs) {
        allr <- allr[1:(nbObs-1)]
        warning(paste("Only first",nbObs-1,"radii used"))
    } else if(length(allr)<(nbObs-1)) {
        stop("'allr' should be of length 1 or nbObs-1")
    }
    
    xy <- simLG_rcpp(nbObs=nbObs, beta=beta, allr=allr, cov=cov, xy0=xy0, 
                     lim=lim, res=res, norm=ifelse(norm,1,0), npts=npts)
    return(xy)
}
