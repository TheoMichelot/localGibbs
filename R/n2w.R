
#' Transform parameters from natural to working scale
#' 
#' @param beta Coefficients of the RSF
#' @param rdist Distribution of the availability radius ("fixed", "exp", or "gamma")
#' @param r Availability radius, if rdist="fixed"
#' @param shape Shape parameter, if rdist="gamma"
#' @param rate Rate parameter, if rdist="exp" or rdist="gamma"
#' @param xy Matrix of observed locations, needed to derive maximum step length
#' if rdist="fixed"
#' 
#' @return Vector of parameters on the working scale
#' 
#' @export
n2w <- function(beta, rdist=c("fixed", "exp", "gamma"), r=NULL, 
                shape=NULL, rate=NULL, xy=NULL)
{
    if(rdist=="fixed") {
        if(is.null(r) | is.null(xy))
            stop("'r' and 'xy' must be provided if rdist='fixed'")
        steps <- sqrt(rowSums((xy[-1,]-xy[-nrow(xy),])^2))
        stepmax <- max(steps, na.rm=TRUE)
        wpar <- c(beta, log(r - stepmax/2))
    } else if(rdist=="exp") {
        if(is.null(rate))
            stop("'rate' must be provided if rdist='exp'")
        wpar <- c(beta, log(rate))
    } else if(rdist=="gamma") {
        if(is.null(shape) | is.null(rate))
            stop("'shape' and 'rate' must be provided if rdist='gamma'")
        wpar <- c(beta, log(shape), log(rate))
    }
    return(wpar)
}
