
#' Fit RSF
#' 
#' @param xy Matrix of observed locations
#' @param method Method to fit the RSF, either "poisson" for a
#' Poisson GLM, or "logit" for a logistic regression
#' @param covlist List of covariate rasters
#' @param nbin Number of bins in each dimension for Poisson GLM
#' @param nzeros Number of zeros to sample for logistic regression
#' 
#' @return Output of glm
#' 
#' @importFrom ash bin2
#' @importFrom reshape2 melt
#' @importFrom stats binomial glm
#' @importFrom raster extent extract
#' @export
fitRSF <- function(xy, method=c("poisson","logit"), covlist, nbin=NULL, nzeros=NULL)
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
        
        modRSF <- glm(countsLG~., data=RSFdata, family="poisson")
    } else if(method=="logit") {
        xyzeros <- matrix(c(runif(nzeros,lim[1],lim[2]), runif(nzeros,lim[3],lim[4])), ncol=2)
        
        # extract cov values in cells
        allcov <- lapply(covlist, extract, y=rbind(xy,xyzeros))
        RSFdata <- as.data.frame(allcov)
        colnames(RSFdata) <- NULL
        RSFdata <- cbind(used=c(rep(1,nrow(xy)), rep(0,nzeros)), RSFdata)
        
        modRSF <- glm(used~., family=binomial(link="logit"), data=RSFdata)
    }
    
    return(modRSF)    
}
