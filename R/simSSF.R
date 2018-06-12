
#' Simulate from SSF
#' 
#' @param nbObs Number of observation
#' @param beta Vector of coefficients of the SSF
#' @param xy1 First location
#' @param xy2 Second location
#' @param nzeros Number of zeros to sample for each step
#' @param covlist List of covariate rasters
#' @param stepPar Shape and rate of gamma distribution
#' @param anglePar Mean and concentration of von Mises distribution
#' 
#' @importFrom CircStats rvm
#' @export
simSSF <- function(nbObs, beta, xy1, xy2, nzeros, covlist, stepPar, anglePar)
{
    cov <- listToArray(covlist)
    lim <- as.vector(extent(covlist[[1]]))
    res <- c(xres(covlist[[1]]),yres(covlist[[1]]))
    
    xy <- matrix(NA, nbObs,2)
    xy[1,] <- xy1
    xy[2,] <- xy2

    t <- 3
    cat("Simulation from SSF...\n")
    pb <- txtProgressBar(min = 1, max = nbObs, style = 3)
    while(t<=nbObs) {
        gridl <- rgamma(nzeros, stepPar[1], stepPar[2])
        gridth <- rvm(nzeros, anglePar[1], anglePar[2])
        gridbear <- gridth + atan2(xy[t-1,2] - xy[t-2,2], xy[t-1,1]-xy[t-2,1])
        grid <- rep(xy[t-1,],each=nzeros) + gridl * cbind(cos(gridbear),sin(gridbear))
        
        allrsf <- rsfvec(grid,beta,cov,lim,res)
        
        if(sum(allrsf)>0) {
            ind <- sample(1:nzeros, size=1, prob=allrsf/sum(allrsf))
            xy[t,] <- grid[ind,]
            t <- t+1
        }
        setTxtProgressBar(pb, t)
    }
    cat("\n")
    
    return(xy)
}
