/*
* LGIv Model
* -----------------------------------------------------
* Copyright: J. A. Formaggio <josephf@mit.edu>
*
* Date: 10 July 2015
*
* Purpose: 
*
*		Post-processing of LGIv violation analysis
*
* Collaboration:  J. A. Formaggio, M. Murskyj, T. Weiss, D. Kaiser
*/

data {

//   Number of data inputs 

     int nData;
     int nPoints[nData];
     int nLGIv_MC[nData];
     int nLGIv[nData];

}

transformed data {

     int nP;
     int nMeasured;
     int nDist[nPoints[1]];
     int nZero;
     int nTotal;     

     nMeasured <- nLGIv[1];
     nP <- nPoints[1];
     
     nZero <- 0;
     for (i in 1:nP) nDist[i] <- 0;
     for (i in 1:nData) {
     	 if (nLGIv_MC[i] > 0) {
	    nDist[nLGIv_MC[i]] <- nDist[nLGIv_MC[i]] + 1;
	 } else {
	    nZero <- nZero + 1;
	 }
     }
     nTotal <- sum(nDist) + nZero;

}

parameters {

    real<lower=-1.0,upper=1.0> pQ;
    real<lower=0.> alpha;
    real<lower=0.> beta;

}


model {
    
//  Fit that distribution to a Gaussian x Binomial distribution.
      
      vector[nP] mu;
      real nNorm;
      
      nNorm <- (1.+pQ) * nTotal;
      for (i in 1:nP) {
          mu[i] <- nNorm * exp(beta_binomial_log(i, nP, alpha, beta));
      }
      
      nDist ~ poisson(mu);
      nZero ~ poisson(nNorm * exp(beta_binomial_log(0, nP, alpha, beta)));
}

generated quantities{

     int n_recon;
     real survival_prob;

     n_recon <- beta_binomial_rng(nP, alpha, beta);
     survival_prob <- beta_binomial_cdf(nMeasured, nP, alpha, beta);

}