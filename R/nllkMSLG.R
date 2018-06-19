
#' Negative log-likelihood for multistate local Gibbs model
#' 
#' @param par Vector of parameters on working scale
#' @param xy Matrix of observed locations
#' @param nstate Number of states
#' @param MCgrids List of Monte Carlo samples, with components gridc and gridz. 
#' As output by \code{\link{MCsample}}.
#' @param cov Array of covariates (one layer for each covariate)
#' @param lim Limits of the covariate rasters.
#' @param res Resolution of the covariate rasters.
#' @param norm Logical. TRUE if normal transition density.
#' 
#' @return Negative log-likelihood
#' 
#' @export
nllkMSLG <- function(par, xy, nstate, MCgrids, cov, lim, res, norm=FALSE)
{
    nobs <- nrow(xy)
    steps <- sqrt(rowSums((xy[-1,]-xy[-nrow(xy),])^2))
    
    # unpack parameters
    ncov <- dim(cov)[3]
    npar <- w2n(par, rdist="multistate", nstate=nstate, xy=xy, norm=norm)
    beta <- npar$beta
    r <- NULL
    sigma <- NULL
    if(!norm)
        r <- npar$r
    else
        sigma <- npar$sigma
    gamma <- npar$gamma
    
    # unpack Monte Carlo samples
    gridc <- MCgrids$gridc
    gridz <- MCgrids$gridz
    
    # calculate state dependent probs
    p <- HMMprobs(xy=xy, r=r, sigma=sigma, beta=beta, gridc=gridc, gridz=gridz, 
                  cov=cov, lim=lim, res=res, norm=norm)
    
    if(any(p[,1]==0 & p[,2]==0)) {
        warning(paste("At least one matrix of state-dependent densities is",
                      "zero. Returning nllk=1e8."))
        return(1e8)
    }
    
    # initial distribution
    delta <- rep(1,nstate)/nstate
    
    # forward algorithm
    llk <- 0
    foo <- delta*p[1,]
    for(t in 2:(nobs-1)) {
        foo <- foo%*%gamma*p[t,]
        sumfoo <- sum(foo)
        llk <- llk + log(sumfoo)
        if(any(is.nan(foo/sumfoo))) {
            cat("foo =",foo,", p[t,] =",p[t,],"\n")
            stop("Error in forward algorithm")
        }
        foo <- foo/sumfoo
    }
    
    return(-llk)
}
