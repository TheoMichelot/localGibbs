
#' Local Gibbs negative log likelihood (normal transition density)
#' 
#' @param par Vector of parameters on working scale
#' @param ID Vector of track ID
#' @param xy Matrix of observed locations
#' @param MCgrids List of Monte Carlo samples for the pseudo-likelihood evaluation
#' @param cov Array of covariates (one layer for each covariate)
#' @param lim Limits of the covariate rasters.
#' @param res Resolution of the covariate rasters.
#' 
#' @export
nllkLG_norm <- function(par, ID=NULL, xy, dt=NULL, MCgrids, cov, lim, res)
{
    # consider unique track if ID==NULL
    if(is.null(ID))
        ID <- rep(1,nrow(xy))
    
    # consider regular time intervals if dt==NULL
    if(is.null(dt))
        dt <- rep(1, nrow(xy)-1)
    
    if(length(par)!=dim(cov)[3]+1)
        stop("Length of 'par' incompatible with dimensions of 'cov'")
    
    # unpack parameters
    beta <- par[1:dim(cov)[3]]
    sigma <- exp(par[dim(cov)[3]+1])
    
    nllk <- nllkLG_norm_rcpp(beta=beta, sigma=sigma, ID=ID, xy=xy, dt=dt, gridc=MCgrids$gridc,
                             gridz=MCgrids$gridz, cov=cov, lim=lim, res=res)$nllk
    
    return(nllk)
}
