#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "scalez.hpp"
#include "rsf.hpp"

//' Negative log-likelihood for the local Gibbs model (C++)
//'
//' @param beta Parameters of the RSF
//' @param sigma Standard deviation parameter
//' @param ID Vector of track IDs
//' @param xy Matrix of observed locations
//' @param dt Vector of time intervals
//' @param gridc Grid for Monte Carlo integration
//' @param gridz Grid for Monte Carlo integration
//' @param cov Array of covariates (one layer for each covariate)
//' @param lim Limits of the covariate rasters.
//' @param res Resolution of the covariate rasters.
// [[Rcpp::export]]
Rcpp::List nllkLG_norm_rcpp(arma::vec beta, double sigma, arma::vec ID, arma::mat xy,
                            arma::vec dt, arma::mat gridc, arma::mat gridz, arma::cube& cov, 
                            arma::vec lim, arma::vec res)
{
    int nbObs = xy.n_rows;
    
    // lengths of steps
    arma::vec steps(nbObs-1);
    for(int i=0; i<nbObs-1; i++)
        steps(i) = sqrt(pow(xy(i+1,0)-xy(i,0),2) + pow(xy(i+1,1)-xy(i,1),2));
    
    int nc = gridc.n_rows;
    int nz = gridz.n_rows;
    
    arma::vec logpall(nbObs-1);
    double llk = 0;
    for(int t=1; t<nbObs; t++) {
        // no contribution if first location of a track or if missing data
        if(ID(t)==ID(t-1) & arma::is_finite(steps(t-1))) {
            
            double scaled_sigma = sqrt(dt(t-1)) * sigma;
            
            // scale samples to N(0,sigma^2)
            arma::mat gridc_sc = gridc * scaled_sigma;
            arma::mat gridz_sc = gridz * scaled_sigma;
            
            arma::mat gridc_tr(nc,2);
            gridc_tr.col(0) = gridc_sc.col(0) + xy(t-1,0);
            gridc_tr.col(1) = gridc_sc.col(1) + xy(t-1,1);
            
            double sumc = 0;
            int count = nc;
            for(int ic=0; ic<nc; ic++) {
                
                arma::mat gridz_tr(nz,2);
                gridz_tr.col(0) = gridz_sc.col(0) + gridc_tr(ic,0);
                gridz_tr.col(1) = gridz_sc.col(1) + gridc_tr(ic,1);
                
                double sumz = 0;
                for(int iz=0; iz<nz; iz++) {
                    sumz = sumz + rsf(gridz_tr.row(iz),beta,cov,lim,res);
                }
                
                if(sumz>0) {
                    sumc = sumc + R::dnorm(xy(t,0), gridc_tr(ic,0), scaled_sigma, 0) * 
                        R::dnorm(xy(t,1), gridc_tr(ic,1), scaled_sigma, 0) / sumz;    
                } else {
                    count = count - 1;
                }
            }
            
            double logp = log(nz) - log(count) + log(rsf(xy.row(t),beta,cov,lim,res)) + log(sumc);
            
            logpall(t-1) = logp;
            llk = llk + logp;
        }
    }
    
    return Rcpp::List::create(Rcpp::Named("nllk") = -llk, 
                              Rcpp::Named("logp") = logpall);
}
