
#' Generate Monte Carlo samples
#' 
#' @param nr Size of Monte Carlo sample over r
#' @param nc Size of Monte Carlo sample over c
#' @param nz Size of Monte Carlo sample over z
#' @param norm Logical. If TRUE, standard normal bivariate
#' samples are returned for c and z.
#' 
#' @return A list of three arrays of Monte Carlo samples:
#' gridr (vector), gridc (matrix), and gridz (matrix).
#' 
#' @importFrom lhs randomLHS
#' @export
MCsample <- function(nr=NULL, nc, nz, norm=FALSE)
{
    if(norm) {
        gridr <- NULL
        gridc <- qnorm(randomLHS(nc,2))
        gridz <- qnorm(randomLHS(nz,2))
    } else {
        if(is.null(nr) || nr==0)
            gridr <- NULL
        else
            gridr <- as.vector(randomLHS(nr,1))
        gridc <- randomLHS(nc,2) - 0.5
        gridzfoo <- randomLHS(nz,2)
        gridz <- cbind(sqrt(gridzfoo[,1]), 2*pi*(gridzfoo[,2])-pi)        
    }
    
    return(list(gridr=gridr, gridc=gridc, gridz=gridz))
}
