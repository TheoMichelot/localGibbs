
#' Viterbi algorithm for multistate local Gibbs model
#' 
#' @param xy Matrix of observed locations
#' @param r Vector of availability radii (if norm=FALSE)
#' @param sigma Vector of standard deviations (if norm=TRUE)
#' @param beta Vector of selection parameters
#' @param gamma Transition probability matrix
#' @param MCgrids List of Monte Carlo samples, with components gridc 
#' and gridz. As output by \code{\link{MCsample}}.
#' @param covlist List of covariate rasters
#' @param norm Logical. TRUE if normal transition density. (Only for multistate case)
#' 
#' @export
viterbi <- function(xy, r=NULL, sigma=NULL, beta, gamma, 
                    MCgrids, covlist, norm=FALSE)
{
    nobs <- nrow(xy) - 1
    if(norm)
        nstate <- length(sigma)
    else
        nstate <- length(r)
    
    # unpack arguments
    gridc <- MCgrids$gridc
    gridz <- MCgrids$gridz
    cov <- listToArray(covlist)
    lim <- as.vector(extent(covlist[[1]]))
    res <- c(xres(covlist[[1]]), yres(covlist[[1]]))
    
    # state-dependent probs
    p <- HMMprobs(xy=xy, r=r, sigma=sigma, beta=beta, gridc=gridc, gridz=gridz, 
                  cov=cov, lim=lim, res=res, norm=norm)
    
    # initial distribution
    delta <- rep(1,nstate)/nstate
    
    xi <- matrix(NA,nobs,nstate)
    foo <- delta*p[1,]
    xi[1,] <- foo/sum(foo)
    for(i in 2:nobs) {
        foo <- apply(xi[i-1,]*gamma,2,max)*p[i,]
        xi[i,] <- foo/sum(foo)
    }
    
    stSeq <- rep(NA,nobs)
    stSeq[nobs] <- which.max(xi[nobs,])
    for(i in (nobs-1):1)
        stSeq[i] <- which.max(gamma[,stSeq[i+1]]*xi[i,])
    
    return(stSeq)
}
