
#' Compare two utilisation distributions
#' 
#' @param beta1 Coefficients of the first UD
#' @param beta2 Coefficients of the second UD
#' @param covlist List of covariate rasters
#' 
#' @export
#' @importFrom graphics plot abline
compareUD <- function(beta1, beta2, covlist)
{
    if(length(beta1)!=length(covlist) | length(beta2)!=length(covlist))
        stop("'beta1', 'beta2' and 'covlist' must be of same length")
    
    r1 <- rastRSF(beta=beta1, covlist=covlist)
    s <- sum(values(r1)) * xres(r1) * yres(r1)
    values(r1) <- values(r1)/s
    r2 <- rastRSF(beta=beta2, covlist=covlist)
    s <- sum(values(r2)) * xres(r2) * yres(r2)
    values(r2) <- values(r2)/s
    
    l1 <- min(values(r1),values(r2))
    l2 <- max(values(r1),values(r2))
    plot(values(r1), values(r2), pch=20, cex=0.3, xlim=c(l1,l2), ylim=c(l1,l2))
    abline(0,1)
}
