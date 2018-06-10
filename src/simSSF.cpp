#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "simZeros.hpp"
#include "rsf.hpp"
#include "sample.hpp"

//' SSF simulation function (Fortin et al 2005)
//' 
//' @param nbObs Number of observations
//' @param beta Vector of resource selection coefficients
//' @param xy1 First location
//' @param xy2 Second location
//' @param nzeros Number of zeros (unused steps) to sample for each step
//' @param cov Array of covariates (one layer for each covariate)
//' @param lim Limits of map
//' @param res Resolution of map
//' @param stepprobs Vector of probabilities for intervals of the range of step lengths
//' @param stepbreaks Vector of bounds for intervals of the range of step lengths
//' @param angleprobs Vector of probabilities for intervals of the range of turning angles
//' @param anglebreaks Vector of probabilities for intervals of the range of turning angles
//' 
//' @return Matrix of simulated locations
//' @export
// [[Rcpp::export]]
arma::mat simSSF_rcpp(int nbObs, arma::vec beta, arma::rowvec xy1, arma::rowvec xy2, int nzeros,
                      arma::cube cov, arma::vec lim, arma::vec res, arma::vec stepprobs, 
                      arma::vec stepbreaks, arma::vec angleprobs, arma::vec anglebreaks)
{
    arma::mat xy(nbObs,2);
    xy.zeros();
    xy.row(0) = xy1;
    xy.row(1) = xy2;
    
    arma::mat grid(nzeros,2);
    arma::vec allrsf(nzeros);
    double d, a;
    
    int t = 2;
    while(t<nbObs) {
        // sample nzeros points from empirical step and turn distributions
        grid = simZeros_1step(nzeros, xy.rows(t-2,t-1), stepprobs, stepbreaks, angleprobs, anglebreaks);
        
        for(int i=0; i<nzeros; i++)
            allrsf(i) = rsf(grid.row(i),beta,cov,lim,res);
        
        if(sum(allrsf)>0) {
            int ind = sample(allrsf/sum(allrsf));
            xy.row(t) = grid.row(ind);
            t = t+1;
        }
    }
    
    return xy;
}
