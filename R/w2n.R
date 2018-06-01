
#' Transform parameters from working to natural scale
#' 
#' @param wpar Vector of parameters on working scale
#' @param rdist Distribution of the availability radius ("fixed", "exp", or "gamma")
#' @param xy Matrix of observed locations, needed to derive maximum step length
#' if rdist="fixed"
#' 
#' @return Vector of parameters on the working scale
#' 
#' @export
w2n <- function(wpar, rdist=c("fixed", "exp", "gamma"), xy=NULL)
{
    if(rdist=="fixed") {
        steps <- sqrt(rowSums((xy[-1,]-xy[-nrow(xy),])^2))
        stepmax <- max(steps, na.rm=TRUE)
        beta <- wpar[1:(length(wpar)-1)]
        r <- stepmax/2 + exp(wpar[length(wpar)])
        shape <- NULL
        rate <- NULL
    } else if(rdist=="exp") {
        beta <- wpar[1:(length(wpar)-1)]
        rate <- exp(wpar[length(wpar)])
        r <- NULL
        shape <- NULL
    } else if(rdist=="gamma") {
        beta <- wpar[1:(length(wpar)-2)]
        shape <- exp(wpar[length(wpar)-1])
        rate <- exp(wpar[length(wpar)])
        r <- NULL
    }
    
    return(list(beta=beta, r=r, shape=shape, rate=rate))
}
