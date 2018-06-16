
#' Transform parameters from working to natural scale
#' 
#' @param wpar Vector of parameters on working scale
#' @param rdist Distribution of the availability radius ("fixed", "multistate", "exp", or "gamma")
#' @param nstate Number of states, if rdist="multistate"
#' @param xy Matrix of observed locations, needed to derive maximum step length
#' if rdist="fixed"
#' 
#' @return Vector of parameters on the working scale
#' 
#' @export
w2n <- function(wpar, rdist=c("fixed", "multistate", "exp", "gamma"), nstate=1, xy=NULL)
{
    if(rdist=="fixed") {
        if(is.null(xy))
            stop("'xy' must be provided if rdist='fixed'")
        
        steps <- sqrt(rowSums((xy[-1,]-xy[-nrow(xy),])^2))
        stepmax <- max(steps, na.rm=TRUE)
        
        beta <- wpar[1:(length(wpar)-1)]
        r <- stepmax/2 + exp(wpar[length(wpar)])
        shape <- NULL
        rate <- NULL
        gamma <- NULL
        
    } else if (rdist=="multistate") {
        if(is.null(xy))
            stop("'xy' must be provided if rdist='multistate'")
        if(nstate==1)
            stop("'nstate' should be >=2")
        
        steps <- sqrt(rowSums((xy[-1,]-xy[-nrow(xy),])^2))
        stepmax <- max(steps, na.rm=TRUE)
        
        ncov <- length(wpar) - nstate - nstate*(nstate-1)
        
        beta <- wpar[1:ncov]
        r <- exp(wpar[(ncov+1):(ncov+nstate)])
        # constrain last radius to be larger than half longest step length
        r[nstate] <- r[nstate] + stepmax/2
        gamma <- diag(nstate)
        gamma[!gamma] <- exp(wpar[(ncov+nstate+1):(ncov+nstate*nstate)])
        gamma <- gamma/rowSums(gamma)
        
        shape <- NULL
        rate <- NULL
        
    } else if(rdist=="exp") {
        beta <- wpar[1:(length(wpar)-1)]
        rate <- exp(wpar[length(wpar)])
        r <- NULL
        shape <- NULL
        gamma <- NULL
        
    } else if(rdist=="gamma") {
        beta <- wpar[1:(length(wpar)-2)]
        shape <- exp(wpar[length(wpar)-1])
        rate <- exp(wpar[length(wpar)])
        r <- NULL
        gamma <- NULL
    }
    
    return(list(beta=beta, r=r, shape=shape, rate=rate, gamma=gamma))
}
