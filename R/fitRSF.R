
#' Fit RSF
#' 
#' @param xy Matrix of observed locations
#' @param method Method to fit the RSF, either "poisson" for a
#' Poisson GLM, or "logit" for a logistic regression
#' @param covlist List of covariate rasters
#' @param nbin Number of bins in each dimension for Poisson GLM
#' @param nzeros Number of zeros to sample for logistic regression
#' @param buffsize Size of buffer around the observed locations from which
#' control locations should be sampled, if method='logit'.
#' 
#' @return List of \code{mod} (output of glm) and \code{xyzeros} (matrix of
#' unused locations)
#' 
#' @importFrom ash bin2
#' @importFrom reshape2 melt
#' @importFrom stats binomial glm
#' @importFrom raster extent extract
#' @export
fitRSF <- function(xy, method=c("poisson","logit"), covlist, nbin=NULL, nzeros=NULL, buffsize=NULL)
{
    lim <- as.vector(extent(covlist[[1]]))
 
    if(method=="poisson") {
        # define grid for Poisson GLM
        midx <- lim[1] + (lim[2]-lim[1])/nbin * (0:(nbin-1)+0.5)
        midy <- lim[3] + (lim[4]-lim[3])/nbin * (0:(nbin-1)+0.5)
        xygrid <- as.matrix(expand.grid(midx,midy))
        
        # extract cov values in cells
        gridcov <- lapply(covlist, extract, y=xygrid)
        RSFdata <- as.data.frame(gridcov)
        colnames(RSFdata) <- NULL
        
        # count points in cells
        binsLG <- bin2(x = xy[!is.na(xy[,1]),], 
                       ab = matrix(lim,2,2,byrow=TRUE), 
                       nbin = c(nbin,nbin))
        countsLG <- melt(binsLG$nc)$value
        RSFdata <- cbind(countsLG=countsLG, RSFdata)
        
        mod <- glm(countsLG~., data=RSFdata, family="poisson")
        xyzeros <- NULL    
    } else if(method=="logit") {
        
        if(is.null(buffsize)) {
            xyzeros <- matrix(c(runif(nzeros,lim[1],lim[2]), runif(nzeros,lim[3],lim[4])), ncol=2)            
        } else {
            xyzeros <- matrix(NA, nzeros, 2)
            obslim <- c(range(xy[,1],na.rm=TRUE), range(xy[,2],na.rm=TRUE))
            k <- 1
            while(k <= nzeros) {
                if(k%%100==0)
                    cat("\rSampling zeros... ", k, "/", nzeros, sep="")
                z <- runif(2, c(obslim[1], obslim[3])-buffsize, c(obslim[2], obslim[4])+buffsize)
                d <- sqrt(colSums((t(xy)-z)^2))
                if(min(d,na.rm=TRUE) < buffsize) {
                    xyzeros[k,] <- z
                    k <- k + 1
                }
            }
            cat("\n")
        }
        
        # extract cov values in cells
        allcov <- lapply(covlist, extract, y=rbind(xy,xyzeros))
        RSFdata <- as.data.frame(allcov)
        colnames(RSFdata) <- NULL
        RSFdata <- cbind(used=c(rep(1,nrow(xy)), rep(0,nzeros)), RSFdata)
        
        mod <- glm(used~., family=binomial(link="logit"), data=RSFdata)
    }
    
    return(list(mod=mod,xyzeros=xyzeros))    
}
