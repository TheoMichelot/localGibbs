
#' State-dependent distributions for multistate local Gibbs
#' 
#' @param xy Matrix of observed locations
#' @param r Vector of availability radii 
#' (one for each state -- if norm=FALSE)
#' @param sigma Vector of standard deviations 
#' (one for each state -- if norm=TRUE)
#' @param beta Vector of habitat selection parameters
#' @param gridc Monte Carlo sample of intermediate centres c
#' @param gridz Monte Carlo sample of endpoints z
#' @param cov Array of covariates (one layer for each covariate)
#' @param lim Limits of the covariate rasters
#' @param res Resolution of the covariate rasters
#' @param norm Logical. If TRUE, uses normal densities
#' 
#' @return Matrix with one row for each observed step, and one column
#' for each state.
#' 
#' @importFrom stats dnorm
#' 
#' @export
HMMprobs <- function(xy, r=NULL, sigma=NULL, beta, gridc, gridz, cov, lim, res, norm=FALSE)
{
    nobs <- nrow(xy)
    if(norm)
        nstate <- length(sigma)
    else 
        nstate <- length(r)
    steps <- sqrt(rowSums((xy[-nrow(xy),]-xy[-1,])^2))
    nc <- nrow(gridc)
    nz <- nrow(gridz)
    
    allprobs <- matrix(1, nobs-1, nstate)

    for(state in 1:nstate) {
        if(!norm) {
            if(is.null(r))
                stop("'r' must be provided if norm=FALSE")
            
            # indices of impossible steps (longer than 2r)
            ind0 <- which(steps >= (2*r[state]))
            # indices of possible steps
            ind1 <- which(steps < (2*r[state]) & !is.na(steps))
            
            # p=0 if impossible step
            allprobs[ind0,state] <- 0
            
            # areas of intersection of discs of radius r and centres x_t and x_{t+1}
            Ainter <- rep(0, nobs-1)
            Ainter[ind1] <- 2*r[state]^2*acos(steps[ind1]/(2*r[state])) - r[state]*steps[ind1] * 
                sqrt(1-steps[ind1]^2/(4*r[state]^2));
            
            # loop over possible steps
            for(t in ind1) {
                # scale sample of endpoints
                scaledz <- scalez(gridc=gridc, gridz=gridz, r=r[state], xy0=xy[t,], xy1=xy[t+1,])
                allrsf <- rsfvec(scaledz, beta, cov, lim, res)
                count <- nrow(scaledz)/nz
                
                # sum over z
                sumz <- sapply(1:count, function(i) sum(allrsf[((i-1)*nz+1):(nz*i)]))
                # sum over c
                keep <- which(sumz>0)
                count <- length(keep)
                sumc <- sum(1/sumz[keep])
                
                # pseudo-likelihood of step
                allprobs[t,state] <- rsf(xy[t+1,], beta, cov, lim, res)/(pi*r[state]^2)^2 * 
                    Ainter[t] * nz/count * sumc
            }
        } else {
            if(is.null(sigma))
                stop("'sigma' must be provided if norm=TRUE")
            
            gridc_sc <- sigma[state]*gridc
            gridz_sc <- sigma[state]*gridz
            
            ind <- which(!is.na(steps))
            
            for(t in ind) {
                # translate sample of intermediate points
                allc <- rep(xy[t,], each=nc) + gridc_sc
                # translate sample of endpoints
                allz <- matrix(rep(xy[t,], each=nc*nz), ncol=2)
                allz <- allz + rep(gridc_sc, each=nz)
                allz[,1] <- allz[,1] + rep(gridz_sc[,1], nc)
                allz[,2] <- allz[,2] + rep(gridz_sc[,2], nc)

                allrsf <- rsfvec(allz, beta, cov, lim, res)
                
                # sum over z
                sumz <- sapply(1:nc, function(i) sum(allrsf[((i-1)*nz+1):(nz*i)]))
                keep <- which(sumz>0)
                count <- length(keep)
                # sum over c
                kern <- dnorm(xy[t+1,1], allc[,1], sigma[state]) * 
                    dnorm(xy[t+1,2], allc[,2], sigma[state])
                sumc <- sum(kern[keep]/sumz[keep])
                
                # pseudo-likelihood of step
                allprobs[t,state] <- rsf(xy[t+1,], beta, cov, lim, res) * nz/count * sumc
            }
        }
    }
    
    return(allprobs)
}
