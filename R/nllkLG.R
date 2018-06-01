
#' Negative log-likelihood function for the local Gibbs model
#'
#' @param par Vector of parameters, i.e. resource selection coefficients and
#' parameters for the distribution of the radius.
#' @param ID Vector of track IDs
#' @param xy Matrix of observed locations
#' @param rdist Distribution of availability radius ("fixed", "exp", or "gamma")
#' @param MCgrids List of Monte Carlo samples, with components gridr (if needed), 
#' gridc, and gridz. As output by \code{\link{MCsample}}.
#' @param cov Array of covariates (one layer for each covariate)
#' @param lim Limits of the covariate rasters.
#' @param res Resolution of the covariate rasters.
#' 
#' @export
#' 
#' @useDynLib localGibbs
nllkLG <- function(par, ID=NULL, xy, rdist=c("fixed", "exp", "gamma"), MCgrids, cov, lim, res)
{
    # consider unique track if ID==NULL
    if(is.null(ID))
        ID <- rep(1,nrow(xy))

    gridc <- MCgrids$gridc
    gridz <- MCgrids$gridz
    
    # unpack parameters
    tpar <- w2n(wpar=par, rdist=rdist, xy=xy)
    beta <- tpar$beta
    r <- tpar$r
    shape <- tpar$shape
    rate <- tpar$rate
    
    if(length(beta)!=dim(cov)[3])
        stop("Wrong number of parameters provided")
    
    if(rdist=="exp") {
        gridr <- MCgrids$gridr
        truncr <- truncgridr(shape, rate, ID, xy, gridr)
    } else if(rdist=="gamma") {
        truncr <- truncgridr(shape, rate, ID, xy, gridr)
    } else if(rdist=="fixed") {
        truncr <- matrix(r, nrow=nrow(xy)-1, ncol=1)
    }
    
    # Rcpp doesn't cope with NULL
    if(is.null(shape))
        shape <- NA
    if(is.null(rate))
        rate <- NA
    
    # compute negative log-likelihood in C++
    nllk <- nllkLG_rcpp(beta=beta, shape=shape, rate=rate, ID=ID, xy=xy, rdist=rdist,
                        truncr=truncr, gridc=gridc, gridz=gridz, cov=cov, lim=lim, res=res)

    return(nllk)
}
