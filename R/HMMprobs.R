
#' State-dependent distributions for multistate local Gibbs
#' 
#' @param xy Matrix of observed locations
#' @param r Vector of availability radii (one for each state)
#' @param beta Vector of habitat selection parameters
#' @param gridc Monte Carlo sample of intermediate centres c
#' @param gridz Monte Carlo sample of endpoints z
#' @param cov Array of covariates (one layer for each covariate)
#' @param lim Limits of the covariate rasters.
#' @param res Resolution of the covariate rasters.
#' 
#' @return Matrix with one row for each observed step, and one column
#' for each state.
#' 
#' @export
HMMprobs <- function(xy, r, beta, gridc, gridz, cov, lim, res)
{
    nobs <- nrow(xy)
    nstate <- length(r)
    steps <- sqrt(rowSums((xy[-nrow(xy),]-xy[-1,])^2))
    nz <- nrow(gridz)
    
    allprobs <- matrix(1, nobs-1, nstate)
    for(state in 1:nstate) {
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
    }
    
    return(allprobs)
}
