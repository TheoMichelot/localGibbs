#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Resource selection function
//'
//' @param xy Location where the RSF should be evaluated
//' @param beta Vector of resource selection coefficients
//' @param cov Array of covariates (one layer for each covariate)
//' @param lim Four values: x min, x max, y min, y max
//' @param res Two values: resolution in x, and resolution in y
// [[Rcpp::export]]
double rsf(arma::rowvec xy, arma::vec beta, arma::cube& cov, arma::vec lim, arma::vec res)
{
    // return 0 if outside region
    if(xy(0)<=lim(0) || xy(0)>=lim(1) || xy(1)<=lim(2)  || xy(1)>=lim(3))
        return 0;
    
    // extract covariate values
    int i = floor((lim(3)-xy(1))/res(1));
    int j = floor((xy(0)-lim(0))/res(0));
    
    arma::vec covxy(size(cov)(2));
    
    if(i>=size(cov)(0) || j>=size(cov)(1)) {
        Rcpp::stop("Error in rsf: index out of bounds. Check lim and res.");
    } else {
        covxy = cov.tube(i,j);
    }
    
    // return 0 if NA covariate values (covariate unmapped)
    for(int k=0; k<covxy.size(); k++) {
        if(!arma::is_finite(covxy(k)))
            return 0;
    }
    
    double reg = sum(covxy%beta);
    return exp(reg);
}

//' Resource selection function (vectorized for xy)
//' 
//' @param xy Matrix of points where the RSF should be evaluated
//' @param beta Vector of resource selection coefficients
//' @param cov Array of covariates (one layer for each covariate)
//' @param lim Four values: x min, x max, y min, y max
//' @param res Two values: resolution in x, and resolution in y
//' 
//' @export
// [[Rcpp::export]]
arma::vec rsfvec(arma::mat xy, arma::vec beta, arma::cube& cov, arma::vec lim, arma::vec res)
{
    int nbObs = xy.n_rows;
    arma::vec allrsf(nbObs);
    
    for(int i=0; i<nbObs; i++)
        allrsf(i) = rsf(xy.row(i), beta, cov, lim, res);
    
    return allrsf;
}
