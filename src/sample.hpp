#ifndef _SAMPLE_
#define _SAMPLE_

// Sample from a vector of probabilities
int sample(arma::vec probs)
{
    double u = R::runif(0,1);
    double q = 0;
    int k = 0;

    while(q<u) {
        q = q + probs(k);
        k = k+1;
    }

    return k-1;
}

#endif
