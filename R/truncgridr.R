
#' Truncate the sample of avilability radii
#'
#' @param rdist Distribution of availability radius r. Can be
#' "exp", "gamma", or "weibull".
#' @param shape Shape parameter of r distribution (NA if r~exp)
#' @param rate Rate parameter of r distribution
#' @param ID Vector of track IDs
#' @param xy Matrix of observed locations
#' @param gridr Grid for Monte Carlo integration
#' 
#' @importFrom truncdist qtrunc
truncgridr <- function(rdist = c("exp","gamma","weibull"), shape=NULL, rate, ID, xy, gridr)
{
    # step lengths
    steps <- sqrt(rowSums((xy[-1,] - xy[-nrow(xy),])^2))

    # no contribution if first obs in track, or if missing data
    allt <- which(ID[-length(ID)]==ID[-1] & !is.na(steps))

    truncr <- matrix(NA,nrow=nrow(xy)-1,ncol=length(gridr))

    arglist <- list(p = gridr,
                    spec = rdist,
                    b = Inf)
    
    if(rdist=="gamma" | rdist=="weibull")
        arglist$shape <- shape
    if(rdist=="exp" | rdist=="gamma")
        arglist$rate <- rate
    if(rdist=="weibull")
        arglist$scale <- rate
    
    for(t in allt) {
        arglist$a <- steps[t]/2

        tryCatch(
            truncr[t,] <- do.call(qtrunc, arglist),
            error = function(e) {
                print(arglist)
                stop(e)
            }
        )
    }

    return(truncr)
}
