
#' Truncate the sample of avilability radii
#'
#' @param shape Shape parameter of r distribution (NA if r~exp)
#' @param rate Rate parameter of r distribution
#' @param ID Vector of track IDs
#' @param xy Matrix of observed locations
#' @param gridr Grid for Monte Carlo integration
#' 
#' @importFrom truncdist qtrunc
truncgridr <- function(shape=NULL, rate, ID, xy, gridr)
{
    # step lengths
    steps <- sqrt(rowSums((xy[-1,] - xy[-nrow(xy),])^2))

    rdist <- ifelse(is.null(shape), "exp", "gamma")

    # no contribution if first obs in track, or if missing data
    allt <- which(ID[-length(ID)]==ID[-1] & !is.na(steps))

    truncr <- matrix(NA,nrow=nrow(xy)-1,ncol=length(gridr))

    for(t in allt) {
        arglist <- list(p = gridr,
                        spec = rdist,
                        a = steps[t]/2,
                        b = Inf,
                        rate = rate)

        if(!is.null(shape))
            arglist$shape <- shape

        tryCatch(
            truncr[t,] <- do.call(qtrunc, arglist),
            error = function(e) stop("Error with qtrunc...")
        )
    }

    return(truncr)
}
