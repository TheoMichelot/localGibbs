
#ifndef _SIMZEROS_
#define _SIMZEROS_

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Generate zeros from step-and-turn empirical distributions, for one observed step
arma::mat simZeros_1step(int nzeros, arma::mat xy, arma::vec stepprobs, arma::vec stepbreaks,
                         arma::vec angleprobs, arma::vec anglebreaks)
{
    arma::mat zeros(nzeros, 2);
    
    int stepint, angleint;
    double stepzero, anglezero;
    double oldbearing, newbearing;
    
    oldbearing = atan2(xy(1,1)-xy(0,1), xy(1,0)-xy(0,0));
    
    for(int i=0; i<nzeros; i++) {
        stepint = sample(stepprobs);
        stepzero = R::runif(stepbreaks(stepint),stepbreaks(stepint+1));
        angleint = sample(angleprobs);
        anglezero = R::runif(anglebreaks(angleint),anglebreaks(angleint+1));
        newbearing = oldbearing + anglezero;
        
        zeros(i,0) = xy(1,0) + stepzero * cos(newbearing);
        zeros(i,1) = xy(1,1) + stepzero * sin(newbearing);
    }
    
    return zeros;
}

//' Sample zeros for SSF
//' 
//' @param nzeros Number of zeros for each observed location
//' @param xy Matrix of observed locations
//' @param stepprobs Vector of probabilities for intervals of the range of step lengths
//' @param stepbreaks Vector of bounds for intervals of the range of step lengths
//' @param angleprobs Vector of probabilities for intervals of the range of turning angles
//' @param anglebreaks Vector of probabilities for intervals of the range of turning angles
// [[Rcpp::export]]
arma::mat simZeros_rcpp(int nzeros, arma::mat xy, arma::vec stepprobs, arma::vec stepbreaks,
                        arma::vec angleprobs, arma::vec anglebreaks)
{
    int nbObs = xy.n_rows;
    arma::mat zeros((nbObs-2)*nzeros,2);
    
    for(int t=2; t<nbObs; t++) {
        zeros.rows((t-2)*nzeros, (t-1)*nzeros-1) = 
            simZeros_1step(nzeros, xy.rows(t-2,t), stepprobs, stepbreaks, angleprobs, anglebreaks);
    }
    
    return zeros;
}

#endif
