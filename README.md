# Local Gibbs model

This R package implements simulation and estimation functions for the local Gibbs movement model described in:

Michelot, T., Blackwell, P.G., Matthiopoulos, J. (in press).  
["Linking resource selection and step selection models for habitat preferences in animals"](https://esajournals.onlinelibrary.wiley.com/doi/abstract/10.1002/ecy.2452).  
_Ecology_. DOI: 10.1002/ecy.2452.

### Example code

The code below simulates a trajectory from the local Gibbs model.

``` R
# install package from Github
devtools::install_github("TheoMichelot/localGibbs")
library(localGibbs)

# simulate three covariate rasters
covlist <- list()
for(i in 1:3)
    covlist[[i]] <- simRaster(rho=20, lim=c(-100,100,-80,80), res=1)

# simulate from local Gibbs model
nbObs <- 1000 # number of locations to simulate
beta <- c(-2,4,2) # RSF coefficients
allr <- 3 # availability radius (affects perception and movement speed)
xy <- simLG(nbObs=nbObs, beta=beta, allr=allr, covlist=covlist)

# plot simulated track on RSF
plotRaster(rastRSF(beta=beta, covlist=covlist), xy=xy)

```