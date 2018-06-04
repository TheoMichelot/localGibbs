
#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "atan2.hpp"
#include "sample.hpp"
#include "simZeros.hpp"

//' Resource selection function
//'
//' @param xy Matrix of points where the RSF should be evaluated
//' @param beta Vector of resource selection coefficients
//' @param cov Array of covariates (one layer for each covariate)
//' @param lim Four values: x min, x max, y min, y max
//' @param res Two values: resolution in x, and resolution in y
// [[Rcpp::export]]
double rsf(arma::rowvec xy, arma::vec beta, arma::cube& cov, arma::vec lim, arma::vec res)
{
    // return 0 if outside region
    if(xy(0)<lim(0) || xy(0)>lim(1) || xy(1)<lim(2)  || xy(1)>lim(3))
        return 0;
    
    // extract covariate values
    int i = floor((lim(3)-xy(1))/res(1));
    int j = floor((xy(0)-lim(0))/res(0));
    arma::vec covxy = cov.tube(i,j);
    
    // return 0 if NA covariate values (covariate unmapped)
    for(int k=0; k<covxy.size(); k++) {
        if(!arma::is_finite(covxy(k)))
            return 0;
    }
    
    double reg = sum(covxy%beta);
    return exp(reg);
}

//' Scale Monte Carlo sample of end points
//'
//' @param gridc (Non-scaled) grid of intermediate centres, as a matrix
//' with two columns.
//' @param gridz (Non-scaled) grid of end points, as a matrix with
//' two columns.
//' @param r Availability radius
//' @param xy0 Origin of step
//' @param xy1 End of step
// [[Rcpp::export]]
arma::mat scalez(arma::mat gridc, arma::mat gridz, double r,
                 arma::rowvec xy0, arma::rowvec xy1)
{
    int nc = gridc.n_rows;
    int nz = gridz.n_rows;
    arma::mat gridc_sc(nc,2);
    
    double step = sqrt(pow(xy1(0)-xy0(0),2) + pow(xy1(1)-xy0(1),2));
    
    // project to rectangle with correct dimensions
    gridc_sc.col(0) = gridc.col(0)*(2*r-step);
    gridc_sc.col(1) = gridc.col(1)*2*sqrt(r*r-step*step/4);
    
    // convert to polar coordinates
    arma::vec l = sqrt(gridc_sc.col(0)%gridc_sc.col(0)+gridc_sc.col(1)%gridc_sc.col(1));
    arma::vec th = atan2vec(gridc_sc.col(1),gridc_sc.col(0));
    
    // rotate
    th = th + atan2(xy1(1)-xy0(1),xy1(0)-xy0(0));
    // convert back to cartesian coordinates
    gridc_sc.col(0) = l % cos(th);
    gridc_sc.col(1) = l % sin(th);
    
    // translation to middle point between x_t-1 and x_t
    gridc_sc.col(0) = gridc_sc.col(0) + (xy0(0)+xy1(0))/2;
    gridc_sc.col(1) = gridc_sc.col(1) + (xy0(1)+xy1(1))/2;
    
    // only keep points from the rectangle which are in the intersection of the discs
    arma::mat subC(nc,2);
    double count = 0;
    for(int j=0; j<nc; j++) {
        bool cond = sqrt(pow(gridc_sc(j,0)-xy0(0),2)+pow(gridc_sc(j,1)-xy0(1),2))<r &&
            sqrt(pow(gridc_sc(j,0)-xy1(0),2)+pow(gridc_sc(j,1)-xy1(1),2))<r;
        
        if(cond) {
            subC.row(count) = gridc_sc.row(j);
            count = count + 1;
        }
    }
    
    // scale gridz to disc of radius r and compute cartesian coordinates
    arma::mat gridz_sc(nz,2);
    for(int j=0; j<nz; j++) {
        gridz_sc(j,0) = r*gridz(j,0) * cos(gridz(j,1));
        gridz_sc(j,1) = r*gridz(j,0) * sin(gridz(j,1));
    }
    
    arma::mat gridz_all(count*nz,2);
    
    for(int ic=0; ic<count; ic++)
        for(int iz=0; iz<nz; iz++)
            gridz_all.row(ic*nz+iz) = subC.row(ic) + gridz_sc.row(iz);
    
    return gridz_all;
}

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
double nllkLG_rcpp(arma::vec beta, double shape, double rate, arma::vec ID, arma::mat xy,
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

            llk = llk + logp;
        }
    }

    return -llk;
}

//' Simulation function (local Gibbs sampler)
//'
//' @param nbObs Number of observations
//' @param beta Vector of resource selection coefficients
//' @param allr Vector of radii for movement kernel (of length nbObs)
//' @param cov Array of covariates (one layer for each covariate)
//' @param xy0 Initial location
//' @param lim Limits of map
//' @param res Resolution of map
//' @export
// [[Rcpp::export]]
arma::mat simLG_rcpp(int nbObs, arma::vec beta, arma::vec allr, arma::cube& cov,
                     arma::rowvec xy0, arma::vec lim, arma::vec res)
{
    if(xy0(0)<lim(0) || xy0(0)>lim(1) || xy0(1)<lim(2) || xy0(1)>lim(3))
        Rcpp::stop("The initial location must be within the limits of the map.");
    
    arma::mat xy(nbObs,2);
    xy.zeros();
    xy.row(0) = xy0;
    
    double d, a;
    arma::rowvec C(2);
    arma::mat grid(100,2);
    arma::vec allrsf(100);
    
    int t = 1;
    while(t<nbObs) {
        d = sqrt(R::runif(0,allr(t)*allr(t)));
        a = R::runif(-M_PI,M_PI);
        // centre of disc
        C(0) = xy(t-1,0) + d*cos(a);
        C(1) = xy(t-1,1) + d*sin(a);
        
        // sample 100 points in the circle
        for(int i=0; i<100; i++) {
            d = sqrt(R::runif(0,allr(t)*allr(t)));
            a = R::runif(-M_PI,M_PI);
            grid(i,0) = C(0) + d * cos(a);
            grid(i,1) = C(1) + d * sin(a);
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
