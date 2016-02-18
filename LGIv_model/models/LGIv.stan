/*
* LGIv Model
* -----------------------------------------------------
* Copyright: J. A. Formaggio <josephf@mit.edu>
*
* Date: 10 July 2015
*
* Purpose: 
*
*		Analyze exisiting data for the purpose of determining LGIv violation
*
* Collaboration:  J. A. Formaggio, M. Murskyj, T. Weiss, D. Kaiser
*/

functions{

// Basic constants

  real c() { return  299792458.;}			   	   	 	// Speed of light in m/s
  real hbar() { return 6.58211928e-16;}			   	   	 	// hbar in eV-s

// Generate the theoretical prediction for K assuming neutrino oscillations
// Uses formula K_n = (n-2) - 2 sin^2(2theta) (\sum_i^{n-1} sin^2(\psi_i,i+1) - sin^2(\psi_n,1))

  real KGenerator(vector phase, real sin_sq_2theta) {

       int n;
       vector[num_elements(phase)] sink;
       real K;

       n <- num_elements(phase);

       for (i in 1:n) sink[i] <- square(sin(phase[i]));       
       K <- (n - 2.0) - 2.0 * sin_sq_2theta * (sum(sink) - 2.0* sink[n]);
       return K;
  }

//  Calculate the kinematic phase for the neutrino oscillation formula

  real get_phase(real deltam2, real x){
       return deltam2 * x / (hbar() * c()) / 4.;
  }

//  Get a classical version of Cij.  Assumes [Q_i, Q_j ] = 0 (Commutation)

  real get_K_classical(int[] phase, vector Pee){

     real K;
     vector[num_elements(phase)-1] Cij;
     
     for (n in 1:num_elements(phase)-1) Cij[n] <- 2.*Pee[phase[n]]-1.0;
     
     K <- sum(Cij) - prod(Cij);
     
     return K;
     
  }

//  Get the quantum-mechanical version of K. Assumes Cij = 2 P(t_i - t_j) - 1

  real get_K(int[] phase, vector Pee){

     real K;
     vector[num_elements(phase)-1] Cij;
     real Cin;

     for (n in 1:num_elements(phase)-1) Cij[n] <- (2.*Pee[phase[n]]-1.0);
     Cin <- (2.*Pee[phase[num_elements(phase)]]-1.0);

     K <- sum(Cij) - Cin;
     
     return K;
     
  }

// Calculate the error on K given n points (uncorrelated, use with caution)

real get_K_err(int[] phase, vector PeeErr){

     real Cijerr;
     real Kerr;

     Kerr <- 0.;
     for (n in 1:num_elements(phase)) {
     	 Cijerr <- 2.*PeeErr[phase[n]];
	 Kerr <- Kerr + Cijerr * Cijerr;
     }

//   Add correlation term if phases are the same

     for (n in 1:num_elements(phase)-1) {
     	 if (phase[n]==phase[n+1]) Kerr <- Kerr + 2.* (2.*PeeErr[phase[n]]) * (2.*PeeErr[phase[n+1]]);
     }
     
     return sqrt(Kerr);
     
  }

//  Given n points find the closest points of the sum of x onto t for a given tolerance t_err

  int find_closest_point(real x, vector t, real t_err){

      int n;
      int iFound;
      vector[num_elements(t)] ldiff ;
      int ldiff_index[num_elements(t)];

      n <- num_elements(t);
      iFound <- -1;
      
      for (i in 1:n) {
      	  ldiff[i] <- fabs(x - t[i])/t_err;
      }
      ldiff_index <- sort_indices_asc(ldiff);
      if (ldiff[ldiff_index[1]] < 1.0) {
      	 iFound <- ldiff_index[1];
      }
      return iFound;

  }

}

data {

//   Order of violation

     int<lower=3,upper=4> nOrder;

//   Choose reality

     int isQuantum;

//   Choose whether to fit oscillation parameters to K or just to oscillation measurement Pee

     int fitK;

//   Number of data inputs 

     int nData;

//   Input L/E, probabilities, and errors

     real Distance;
     real Conversion;
     real relEnergyErr;
     vector[nData] Prob;
     vector[nData] Prob_err;
     vector[nData] Energy;

}


transformed data {

    int nPoints;
    
    vector[nData] L_E;
    vector[nData*nData] K_meas;
    vector[nData*nData] K_meas_sigma;
    
    real K_lower;
    real KC_upper;
    real KQ_upper;
    real K_upper;
    real L_E_tol;
    int phase_array[nData * nData, nOrder];

    int q;

    K_lower <-  2.0-nOrder;
    if (nOrder%2==1) K_lower <- -1.0*nOrder;

    KC_upper <- (nOrder-2.0);
    KQ_upper <- nOrder * cos(pi()/nOrder);
    K_upper <- KC_upper;
    if (isQuantum==1) K_upper <- KQ_upper;

//  Convert from energy to L/E

    L_E <- (Distance * Conversion) ./ Energy;

//  Find the set of data for pairs for either K_3 or K_4

    nPoints <- 0;
    	for (i in 1:nData) {
    	    for (j in i:nData) {
	    	if (nOrder > 3) {
	       	   for (k in j:nData) {
		       L_E_tol <- (L_E[i] + L_E[j] + L_E[k]) * relEnergyErr;
		       q <- find_closest_point(L_E[i] + L_E[j] + L_E[k], L_E, L_E_tol );  
		       if (q > 0) {
		       	  nPoints <- nPoints + 1;
		      	  phase_array[nPoints,1] <- i;  phase_array[nPoints,2] <- j;  phase_array[nPoints,3] <- k;  phase_array[nPoints,4] <- q;
		      	  K_meas[nPoints] <- get_K(phase_array[nPoints],Prob);
		       	  K_meas_sigma[nPoints] <-get_K_err(phase_array[nPoints],Prob_err);
		       	  print(nPoints," : ",Energy[q],"   --   ",K_meas[nPoints] ," +/- ", K_meas_sigma[nPoints]);
		       }
	           }
	        } else {
		   L_E_tol <- (L_E[i] + L_E[j]) * relEnergyErr;
		   q <- find_closest_point(L_E[i] + L_E[j],  L_E, L_E_tol);   
		   if (q > 0) {
		       nPoints <- nPoints + 1;
		       phase_array[nPoints,1] <- i;  phase_array[nPoints,2] <- j;  phase_array[nPoints,3] <- q;
		       K_meas[nPoints] <- get_K(phase_array[nPoints],Prob);
		       K_meas_sigma[nPoints] <-get_K_err(phase_array[nPoints],Prob_err);
		       print(nPoints," : ",Energy[i]," + ", Energy[j]," = ",Energy[i]*Energy[j]/(Energy[i]+Energy[j])," :: ",Energy[q],"   --   ",Prob[i],":",Prob[j],":",Prob[q]," == ",K_meas[nPoints] ," +/- ", K_meas_sigma[nPoints]);
		   }
	        }
	     }
         }

//  Print out the pairs and their measured K

    print("Found ",nPoints," points with limits < ",K_lower,",",K_upper," >");

}

parameters {

    real<lower=-5,upper=0> log10_deltam2;
    real<lower=0.0,upper=1.0> sinsq2theta;
    vector[nData] z1;
    real z2;

    real<lower=0.,upper=1.> pQuantum;
    real<lower=0, upper=1.> pQsigma;
    real<lower=0.,upper=1.> pQ;

}

transformed parameters{

    real deltam2;
    vector[nOrder] phase;
    vector[nData] Pee;
    vector[nData] Pee_sample;
    
    vector[nPoints] K_data;
    vector[nPoints] K_error;
    vector[nPoints] K_sample;
    vector[nPoints] K_MC;

//  Explore delta m^2 in log space

    deltam2 <- pow(10.,log10_deltam2);

//  Calculate the probability of oscillation (Pee or Pmumu)

    for (i in 1:nData){
	Pee[i] <- 1. - sinsq2theta * square(sin(get_phase(deltam2, L_E[i])));
    }

//  Sample the points according ot their measured error.
//  This is ensure proper correlations for different K points is maintained.

    Pee_sample <- Prob + (z1 .* Prob_err);

//  Extract the data for K and its associated (uncorrelated) error

    K_data <- head(K_meas,nPoints);
    K_error <- head(K_meas_sigma, nPoints);

//  Calculate the theoretical prediction of K (K_MC) and the sampled prediction of K (K_sample)
//  Do this for the classical model and the quantum model

    for (i in 1:nPoints) {
	if (isQuantum==1) {
    	   K_sample[i] <- get_K(phase_array[i],Pee_sample);
	   K_MC[i] <- get_K(phase_array[i],Pee);
	} else {
    	   K_sample[i] <- get_K_classical(phase_array[i],Pee_sample);
	   K_MC[i] <- get_K_classical(phase_array[i],Pee);
	}
    }

}

model {

    int nLGIv;

//  Extract normal distributions for all sampled points to construct K

    z1 ~ normal(0.0, 1.0);

//  Chose whether you fit delta_m2 and sinsq2theta based on the survival probability (Pee) or K itself

    if (fitK == 0) {
       Prob ~ normal(Pee, Prob_err);
    } else {
       K_data ~ normal(K_MC, K_error);
    }

//  Determine the number of points that are expected to violate the classical limit.

    nLGIv <- 0;
    for (i in 1:nPoints) nLGIv <- nLGIv + (K_data[i] > KC_upper);
    
//  Fit that distribution to a Gaussian x Binomial distribution.
//  Note:  Because this is slow/unstable, sometimes I just look at the
//         measured violation rather than the sampled variation.

    pQ ~ normal(pQuantum, pQsigma);
    nLGIv ~ binomial(nPoints, pQ);

}

generated quantities {

    int isKfit;
    int iQuantum;
    int nP;
    int nLGIv;
    int nLGIv_MC;
    
    real maxK;
    real minK;
    real deltamuK;
    real sigmaK;

    real maxK_MC;
    real minK_MC;
    real deltamuK_MC;
    real sigmaK_MC;

    real rPick;
    int  nPick;

    real myK_energy[nOrder];
    real myK_phase[nOrder];
    real myK_data;
    real myK_sample;
    real myK_MC;
    real myK_error;
    

//  Output important fit information

    minK <- min(K_sample);
    maxK <- max(K_sample);
    deltamuK  <- mean(K_sample-K_data);
    sigmaK <- sd(K_sample-K_data);
    
    minK_MC <- min(K_MC);
    maxK_MC <- max(K_MC);
    deltamuK_MC  <- mean(K_MC-K_data);
    sigmaK_MC <- sd(K_MC-K_data);

//  Number of observed K points that violate the classical limit.

    nLGIv <- 0;
    for (i in 1:nPoints) {if  (K_data[i] > KC_upper) nLGIv<-nLGIv+1;}	

//  Number of sampled K points that violate the classical limit (the expected violation).

    nLGIv_MC <- 0;
    for (i in 1:nPoints) {
    	if  (K_sample[i] > KC_upper) nLGIv_MC<-nLGIv_MC+1;
    }	

    nP <- nPoints;
    iQuantum <- isQuantum;
    isKfit <- fitK;

//  Now pick a random data point in the sample to make plots, etc.
//  This is done in a bit of a silly way because there is no int<->real conversion in STAN

    rPick <- uniform_rng(1.0,nPoints+1.0);
    for (n in 1:nPoints){
    	if (n >= floor(rPick) && n < ceil(rPick)) nPick <- n;
    }
    
    for (n in 1:nOrder) myK_energy[n] <- Energy[phase_array[nPick,n]];
    myK_data <- K_data[nPick];
    myK_sample <- K_sample[nPick];
    myK_MC <- K_MC[nPick];
    myK_error <- K_error[nPick];
    
}