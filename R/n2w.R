
#' Transform parameters from natural to working scale
#' 
#' @param beta Coefficients of the RSF
#' @param rdist Distribution of the availability radius ("fixed", "multistate", 
#' "exp", or "gamma")
#' @param r Availability radius, if rdist="fixed" (or vector of radii if multistate)
#' @param sigma Vector of standard deviations (if norm=TRUE)
#' @param shape Shape parameter, if rdist="gamma"
#' @param rate Rate parameter, if rdist="exp" or rdist="gamma"
#' @param gamma Transition probability matrix, if rdist="multistate"
#' @param xy Matrix of observed locations, needed to derive maximum step length
#' if rdist="fixed"
#' @param norm Logical. TRUE if normal transition density. (Only for multistate case)
#' 
#' @return Vector of parameters on the working scale
#' 
#' @export
n2w <- function(beta, rdist=c("fixed", "multistate", "exp", "gamma"), r=NULL, sigma=NULL,
                shape=NULL, rate=NULL, gamma=NULL, xy=NULL, norm=FALSE)
{
    if(rdist=="fixed") {
        if(is.null(r) | is.null(xy))
            stop("'r' and 'xy' must be provided if rdist='fixed'")
        steps <- sqrt(rowSums((xy[-1,]-xy[-nrow(xy),])^2))
        stepmax <- max(steps, na.rm=TRUE)
        wpar <- c(beta, log(r - stepmax/2))
        
    } else if(rdist=="multistate") {
        if(is.null(xy) | is.null(gamma))
            stop("'gamma', and 'xy' need to be provided if rdist='multistate'")
        steps <- sqrt(rowSums((xy[-1,]-xy[-nrow(xy),])^2))
        stepmax <- max(steps, na.rm=TRUE)
        
        if(!norm) {
            if(is.null(r))
                stop("'r' must be provided if norm=FALSE")
            
            nstate <- length(r)
            r[length(r)] <- r[length(r)] - stepmax/2
            wr <- log(r)
            foo <- log(gamma/diag(gamma))
            wgamma <- as.vector(foo[!diag(nstate)])
            
            wpar <- c(beta,wr,wgamma)
        } else {
            if(is.null(sigma))
                stop("'sigma' must be provided if norm=TRUE")
            
            nstate <- length(sigma)
            wsigma <- log(sigma)
            foo <- log(gamma/diag(gamma))
            wgamma <- as.vector(foo[!diag(nstate)])
            
            wpar <- c(beta,wsigma,wgamma)
        }
        
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
