#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "rsf.hpp"
#include "sample.hpp"

//' Simulation function (local Gibbs sampler)
//'
//' @param nbObs Number of observations
//' @param beta Vector of resource selection coefficients
//' @param allr Vector of radii for movement kernel (of length nbObs-1)
//' @param cov Array of covariates (one layer for each covariate)
//' @param xy0 Initial location
//' @param lim Limits of map
//' @param res Resolution of map
//' @param norm Logical (0 or 1). If TRUE, a normal transition density is used.
//' @param npts Number of potential endpoints to sample at each time step
//' @export
// [[Rcpp::export]]
arma::mat simLG_rcpp(int nbObs, arma::vec beta, arma::vec allr, arma::cube& cov,
                     arma::rowvec xy0, arma::vec lim, arma::vec res, int norm,
                     int npts)
{
    if(xy0(0)<lim(0) || xy0(0)>lim(1) || xy0(1)<lim(2) || xy0(1)>lim(3))
        Rcpp::stop("The initial location must be within the limits of the map.");
    
    arma::mat xy(nbObs,2);
    xy.zeros();
    xy.row(0) = xy0;
    
    double d, a;
    arma::rowvec C(2);
    arma::mat grid(npts,2);
    arma::vec allrsf(npts);
    
    int t = 1;
    while(t<nbObs) {
        if(norm) {
            C(0) = R::rnorm(xy(t-1,0), allr(t-1));
            C(1) = R::rnorm(xy(t-1,1), allr(t-1));            
        } else {
            d = sqrt(R::runif(0,allr(t-1)*allr(t-1)));
            a = R::runif(-M_PI,M_PI);
            // centre of disc
            C(0) = xy(t-1,0) + d*cos(a);
            C(1) = xy(t-1,1) + d*sin(a);
        }
        
        // sample npts points in the circle
        for(int i=0; i<npts; i++) {
            if(norm) {
                grid(i,0) = R::rnorm(C(0), allr(t-1));
                grid(i,1) = R::rnorm(C(1), allr(t-1));                
            } else {
                d = sqrt(R::runif(0,allr(t-1)*allr(t-1)));
                a = R::runif(-M_PI,M_PI);
                grid(i,0) = C(0) + d * cos(a);
                grid(i,1) = C(1) + d * sin(a);
            }
            allrsf(i) = rsf(grid.row(i),beta,cov,lim,res);
        }
        
        if(sum(allrsf)>0) {
            int ind = sample(allrsf/sum(allrsf));
            xy.row(t) = grid.row(ind);
            t = t+1;
        }
    }
    
    return xy;
}
