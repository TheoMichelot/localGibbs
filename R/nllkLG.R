
#' Negative log-likelihood function for the local Gibbs model
#'
#' @param par Vector of parameters, i.e. resource selection coefficients and
#' parameters for the distribution of the radius.
#' @param ID Vector of track IDs
#' @param xy Matrix of observed locations
#' @param rdist Distribution of availability radius 
#' ("fixed", "exp", "gamma", or "weibull")
#' @param MCgrids List of Monte Carlo samples, with components gridr (if needed), 
#' gridc, and gridz. As output by \code{\link{MCsample}}.
#' @param cov Array of covariates (one layer for each covariate)
#' @param lim Limits of the covariate rasters.
#' @param res Resolution of the covariate rasters.
#' @param debug Logical. If TRUE, the function also returns the log-likelihood
#' of each step
#' 
#' @export
#' 
#' @useDynLib localGibbs
nllkLG <- function(par, ID=NULL, xy, rdist=c("fixed", "exp", "gamma", "weibull"), 
                   MCgrids, cov, lim, res, debug=FALSE)
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
    
    if(rdist=="fixed") {
        truncr <- matrix(r, nrow=nrow(xy)-1, ncol=1)
    } else {
        gridr <- MCgrids$gridr
        truncr <- truncgridr(rdist=rdist, shape=shape, rate=rate, 
                             ID=ID, xy=xy, gridr=gridr)
    }
    
    # Rcpp doesn't cope with NULL
    if(is.null(shape))
        shape <- NA
    if(is.null(rate))
        rate <- NA
    
    # compute negative log-likelihood in C++
    nllk <- nllkLG_rcpp(beta=beta, shape=shape, rate=rate, ID=ID, xy=xy, rdist=rdist,
                        truncr=truncr, gridc=gridc, gridz=gridz, cov=cov, lim=lim, res=res)
    
    if(!debug)
        nllk <- nllk$nllk
    return(nllk)
}
