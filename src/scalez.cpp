#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "atan2.hpp"

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
