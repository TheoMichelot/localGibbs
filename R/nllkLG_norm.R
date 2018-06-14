
#' Local Gibbs negative log likelihood (normal transition density)
#' 
#' @export
nllkLG_norm <- function(par, ID=NULL, xy, MCgrids, cov, lim, res)
{
    # consider unique track if ID==NULL
    if(is.null(ID))
        ID <- rep(1,nrow(xy))
    
    # unpack parameters
    beta <- par[1:dim(cov)[3]]
    sigma <- exp(par[dim(cov)[3]+1])
    
    gridc <- MCgrids$gridc * sigma
    gridz <- MCgrids$gridz * sigma
    
    nllk <- nllkLG_norm_rcpp(beta=beta, sigma=sigma, ID=ID, xy=xy, gridc=gridc, 
                             gridz=gridz, cov=cov, lim=lim, res=res)$nllk    
    
    return(nllk)
}
