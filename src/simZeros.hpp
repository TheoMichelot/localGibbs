#ifndef _SIMZEROS_
#define _SIMZEROS_

arma::mat simZeros_1step(int nzeros, arma::mat xy, arma::vec stepprobs, arma::vec stepbreaks,
                         arma::vec angleprobs, arma::vec anglebreaks);
arma::mat simZeros_rcpp(int nzeros, arma::mat xy, arma::vec stepprobs, arma::vec stepbreaks,
                        arma::vec angleprobs, arma::vec anglebreaks);
    
#endif