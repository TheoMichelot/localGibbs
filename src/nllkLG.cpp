#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "scalez.hpp"
#include "rsf.hpp"

//' Negative log-likelihood for the local Gibbs model (C++)
//'
//' @param beta Parameters of the RSF
//' @param shape Shape parameter of the distribution of r
//' @param rate Rate parameter of the distribution of r
//' @param ID Vector of track IDs
//' @param xy Matrix of observed locations
//' @param rdist Distribution of availability radius ("fixed", "exp", "gamma")
//' @param truncr Truncated grid for Monte Carlo integration
//' @param gridc Grid for Monte Carlo integration
//' @param gridz Grid for Monte Carlo integration
//' @param cov Array of covariates (one layer for each covariate)
//' @param lim Limits of the covariate rasters.
//' @param res Resolution of the covariate rasters.
// [[Rcpp::export]]
Rcpp::List nllkLG_rcpp(arma::vec beta, double shape, double rate, arma::vec ID, arma::mat xy,
                   std::string rdist, arma::mat truncr, arma::mat gridc, arma::mat gridz,
                   arma::cube& cov, arma::vec lim, arma::vec res)
{
    int nbObs = xy.n_rows;
    
    // lengths of steps
    arma::vec steps(nbObs-1);
    for(int i=0; i<nbObs-1; i++)
        steps(i) = sqrt(pow(xy(i+1,0)-xy(i,0),2) + pow(xy(i+1,1)-xy(i,1),2));
    
    int nr = truncr.n_cols;
    int nc = gridc.n_rows;
    int count = 0;
    int nz = gridz.n_rows;
    
    arma::vec logpall(nbObs-1);
    double llk = 0;
    for(int t=1; t<nbObs; t++) {
        // no contribution if first location of a track or if missing data
        if(ID(t)==ID(t-1) & arma::is_finite(steps(t-1))) {
            arma::vec allp(nr);
            
            // loop over radii (and nr=1 if r is fixed)
            for(int i=0; i<nr; i++) {
                double r = truncr(t-1,i);
                
                // obtain grid of z
                arma::mat gridz_all = scalez(gridc, gridz, r, xy.row(t-1), xy.row(t));
                count = gridz_all.n_rows/nz;
                
                // loop over discs
                double sumc = 0;
                for(int ic=0; ic<count; ic++) {
                    double sumz = 0;
                    
                    // loop over points in disc
                    for(int iz=0; iz<nz; iz++)
                        sumz = sumz + rsf(gridz_all.row(ic*nz+iz),beta,cov,lim,res);
                    
                    if(sumz>0)
                        sumc = sumc + 1/sumz;
                }
                
                // area of intersection of discs D_r(x_t) and D_r(x_t-1) (for MC integration)
                double Ainter = 2*r*r*acos(steps(t-1)/(2*r)) - r*steps(t-1) *
                    sqrt(1-steps(t-1)*steps(t-1)/(4*r*r));
                
                allp(i) = Ainter/(count*pow(M_PI,2)*pow(r,4)) * sumc;
            }
            
            double rlogCDF;
            if(rdist=="exp")
                rlogCDF = R::pexp(steps(t-1)/2, 1/rate, 0, 1);
            else if(rdist=="gamma")
                rlogCDF = R::pgamma(steps(t-1)/2, shape, 1/rate, 0, 1);
            else
                rlogCDF = 0;
            
            double logp = log(nz) + log(rsf(xy.row(t),beta,cov,lim,res)) -
                log(nr) + rlogCDF + log(sum(allp));
            
            logpall(t-1) = logp;
            llk = llk + logp;
        }
    }
    
    return Rcpp::List::create(Rcpp::Named("nllk") = -llk, 
                              Rcpp::Named("logp") = logpall);
}
