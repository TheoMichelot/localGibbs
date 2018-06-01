
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
    beta <- par[1:dim(cov)[3]]
    rpar <- exp(par[-(1:dim(cov)[3])])
    
    if(rdist=="exp") {
        if(length(rpar)!=1)
            stop("'par' should be of length ",dim(cov)[3] + 1)
        
        shape <- NA
        rate <- rpar[1]
        rfix <- NA

        gridr <- MCgrids$gridr
        truncr <- truncgridr(shape, rate, ID, xy, gridr)
    } else if(rdist=="gamma") {
        if(length(rpar)!=2)
            stop("'par' should be of length ",dim(cov)[3] + 2)
        
        shape <- rpar[1]
        rate <- rpar[2]
        rfix <- NA

        truncr <- truncgridr(shape, rate, ID, xy, gridr)
    } else if(rdist=="fixed") {
        if(length(rpar)!=1)
            stop("'par' should be of length ",dim(cov)[3] + 1)
        
        shape <- NA
        rate <- NA
        steps <- sqrt(rowSums((xy[-nrow(xy),]-xy[-1,])^2))
        rfix <- rpar[1] + max(steps,na.rm=TRUE)/2

        truncr <- matrix(rfix, nrow=nrow(xy)-1, ncol=1)
    }

    nllk <- nllkLG_rcpp(beta=beta, shape=shape, rate=rate, ID=ID, xy=xy, rdist=rdist,
                        truncr=truncr, gridc=gridc, gridz=gridz, cov=cov, lim=lim, res=res)

    return(nllk)
}
