
#' Transform parameters from working to natural scale
#' 
#' @param wpar Vector of parameters on working scale
#' @param rdist Distribution of the availability radius ("fixed", "multistate", 
#' "exp", "gamma", or "weibull)
#' @param nstate Number of states, if rdist="multistate"
#' @param xy Matrix of observed locations, needed to derive maximum step length
#' if rdist="fixed"
#' @param norm Logical. TRUE if normal transition density. (Only for multistate case)
#' 
#' @return Vector of parameters on the working scale
#' 
#' @export
w2n <- function(wpar, rdist=c("fixed", "multistate", "exp", "gamma", "weibull"), 
                nstate=1, xy=NULL, norm=FALSE)
{
    r <- NULL
    sigma <- NULL
    shape <- NULL
    rate <- NULL
    gamma <- NULL
    
    if(rdist=="fixed") {
        if(is.null(xy))
            stop("'xy' must be provided if rdist='fixed'")
        
        steps <- sqrt(rowSums((xy[-1,]-xy[-nrow(xy),])^2))
        stepmax <- max(steps, na.rm=TRUE)
        
        beta <- wpar[1:(length(wpar)-1)]
        r <- stepmax/2 + exp(wpar[length(wpar)])
        
    } else if (rdist=="multistate") {
        if(is.null(xy))
            stop("'xy' must be provided if rdist='multistate'")
        if(nstate==1)
            stop("'nstate' should be >=2")
        
        steps <- sqrt(rowSums((xy[-1,]-xy[-nrow(xy),])^2))
        stepmax <- max(steps, na.rm=TRUE)
        
        ncov <- length(wpar) - nstate - nstate*(nstate-1)
        beta <- wpar[1:ncov]
        
        if(!norm) {
            r <- exp(wpar[(ncov+1):(ncov+nstate)])
            # constrain last radius to be larger than half longest step length
            r[nstate] <- r[nstate] + stepmax/2
        } else {
            sigma <- exp(wpar[(ncov+1):(ncov+nstate)])
        }

        gamma <- diag(nstate)
        gamma[!gamma] <- exp(wpar[(ncov+nstate+1):(ncov+nstate*nstate)])
        gamma <- gamma/rowSums(gamma)
        
    } else if(rdist=="exp") {
        beta <- wpar[1:(length(wpar)-1)]
        rate <- exp(wpar[length(wpar)])
        
    } else if(rdist=="gamma" || rdist=="weibull") {
        beta <- wpar[1:(length(wpar)-2)]
        shape <- exp(wpar[length(wpar)-1])
        rate <- exp(wpar[length(wpar)])
    }
    
    return(list(beta=beta, r=r, sigma=sigma, shape=shape, rate=rate, gamma=gamma))
}
