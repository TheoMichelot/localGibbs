
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
nllkLG_norm <- function(par, ID=NULL, xy, MCgrids, cov, lim, res)
{
    # consider unique track if ID==NULL
    if(is.null(ID))
        ID <- rep(1,nrow(xy))
    
    if(length(par)!=dim(cov)[3]+1)
        stop("Length of 'par' incompatible with dimensions of 'cov'")
    
    # unpack parameters
    beta <- par[1:dim(cov)[3]]
    sigma <- exp(par[dim(cov)[3]+1])
    
    # scale samples to N(0,sigma^2)
    gridc_sc <- MCgrids$gridc * sigma
    gridz_sc <- MCgrids$gridz * sigma
    
    nllk <- nllkLG_norm_rcpp(beta=beta, sigma=sigma, ID=ID, xy=xy, gridc=gridc_sc,
                             gridz=gridz_sc, cov=cov, lim=lim, res=res)$nllk
    
    return(nllk)
}
