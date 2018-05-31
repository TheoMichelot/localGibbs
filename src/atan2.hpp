#ifndef _ATAN2_
#define _ATAN2_

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// arctangent function with two arguments
double atan2(double x, double y)
{
    if(x>0)
        return atan(y/x);
    else if(y>=0)
        return atan(y/x) + M_PI;
    else
        return atan(y/x) - M_PI;
}

// arctangent function with two arguments (vectorized)
arma::vec atan2vec(arma::vec x, arma::vec y)
{
    arma::vec res(x.size());
    for(int i=0; i<x.size(); i++)
        res(i) = atan2(x(i),y(i));

    return res;
}

#endif
