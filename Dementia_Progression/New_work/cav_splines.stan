data {
  int<lower=0> N; // number of states
  int<lower=0> R; // number of rates modelled 
  int<lower=1> T; // length of data set 
  int<lower=1> K; // number of splines 
  int ID[T]; // id's 
  real Ts[T]; // observed times
  int y[T]; // observations 
  int sex[T]; // Covariates
  real years[T]; // years since transplant
  int years_i[T]; // Integer years since transplant
  matrix[T, K] B1; // Basis at each observation
} 

parameters {
  vector<lower=0, upper=1>[N] p; // misclassification probabilitites
  simplex[N-1] d; // initial distribution 
  vector[R] b_0; //Coefficients
  vector[R] b_1;
  vector[R] b_2;
  vector[K-1] alpha[3]; // spline coefficients 
  vector<lower=0>[3] sigma_a; // coefficient priors sd
}

model {
  vector[N] logp; // log marginal likelihood
  matrix[N, N] Q; // Transition rates matrix, changes for each obs
  matrix[N, N] Q_t; // Auxiliary matrix for Q transposed
  vector[N] logptemp; // auxiliary variable to compute the likelihood
  matrix[N,N] gamma[T]; // Transition probability matrix between t-th obs and previous obs
  matrix[N,N] log_gamma[T]; // log of the ptm
  matrix[N,N] log_gamma_tr[T]; // transposed 
  vector[R] q; // transition rates 
  int per; // auxiliary variable for the number of integer years in each transition
  vector[R] b_1_years[21]; // array with 21 vectors, each to contain b2 times the corresponding integer year
  matrix[N,N] gamma_temp; // auxiliary matrix (computed each integer year) to build ptm 
  vector[N] I; // auxiliary vector to build identity matrix
  matrix[N,N] P; // auxiliary matrix to define the emission distributions
  row_vector[N] nu; // probabilities for first observation
  vector[N] delta_nu; // log probabilities of modified initial distribution
  vector[N] delta; // auxiliary for initial distribution 
  vector[K] asmat[3]; // constrained weights
  
  // define weights
  for(m in 1:3){
    asmat[m] = softmax(alpha[m]);
  }

  // priors for betas
  b_0 ~ normal(-2, 1.1);
  b_1 ~ normal(0, 0.11);
  b_2 ~ normal(0, 0.6);
  
  sigma_a ~ normal(0, 0.5);
  
  // Second-order random walk priors
  for(i in 1:3) {
    for(j in 1:3){
      alpha[j,i] ~ normal(0, 0.5);
    }
  }
  
  for(m in 1:3) {
    alpha[m] ~ normal(0, 0.5);
    for(k in 4:K){
      alpha[m,k] ~ normal(3*alpha[m, k-1]-3*alpha[m,k-2]+alpha[m,k-3], sigma_a[m]);
      }
    }
  
  // resto of the priors
  d[1] ~ normal(0.966, 0.1);
  d[2] ~ normal(0.03, 0.1);
  d[3] ~ normal(0.003, 0.1);
  
  p ~ normal(0.05, .1);

  // populate integer year effect array
  for(i in 1:21){
    b_1_years[i] = (i-8) * b_1;
  } 
  
  // give values to I
  for(i in 1:N){
    I[i] = 1;
  }
  
  // Define emission distribution parameters
  P[1,1] = 1 - p[1];
  P[1,2] = p[1];
  P[1,3] = 0 ;
  P[1,4] = 0 ;
  
  P[2,1] = p[2];
  P[2,2] = 1 - p[2] - p[3];
  P[2,3] = p[3];
  P[2,4] = 0 ;
  
  P[3,1] = 0;
  P[3,2] = p[4];
  P[3,3] = 1 - p[4] ;
  P[3,4] = 0 ;
  
  P[4,1] = 0;
  P[4,2] = 0;
  P[4,3] = 0;
  P[4,4] = 1 ;
  
  // only alive states possible at baseline
  for(i in 1:N-1){ 
    delta[i] = d[i];
  }
  delta[4] = 0;
  
  // derive array of (log-)transition probabilities
  for(t in 1:T) {
    // build NxN identity matrix
    gamma[t] = diag_matrix(I);
    // for first observation
    if( t == 1 || ID[t] != ID[t-1]) {
      // number of integer periods to consider between first obs time and baseline (50)
      per = years_i[t] + 1;
      for(i in 1:per){
        // transitions rates at each year in p 
        q = exp(b_0 + sex[t]*b_2 + b_1_years[i]);
        // transition rate matrix, changes every time 
        Q[1,2] = q[1];
        Q[1,3] = 0;
        Q[1,4] = q[2];
        Q[2,1] = 0;
        Q[2,3] = q[3];
        Q[2,4] = q[4];
        Q[3,1] = 0;
        Q[3,2] = 0;
        Q[3,4] = q[5];
        Q[4,1] = 0;
        Q[4,2] = 0;
        Q[4,3] = 0;
        
  
        Q[1,1] = -sum(q[1:2]);
        Q[2,2] = -sum(q[3:4]);
        Q[3,3] = -q[5];
        Q[4,4] = 0;
        
        // first period 
        if(i == 1){
          if(years[t]<1){ // first observation within the first integer year
            gamma_temp = matrix_exp((years[t])*Q);
          } else{ // integer years within the inter-observations period
            gamma_temp = matrix_exp(Q);
          }
        } else{ // last inger year of the period
          if( i == per){
            gamma_temp = matrix_exp((years[t] - years_i[t])*Q);
          } else{
            gamma_temp = matrix_exp(Q) ;
          }
        } // compute the tpm for all the inter-observations time
        gamma[t] = gamma[t] * gamma_temp ;
      }
    } else{ // when it's not the first observation of the individual
        // number of integer periods to consider between current and previous observation
        per = years_i[t] - years_i[t-1] + 1;
        for(i in 1:per){ // transitions rates at each year in p 
          q = exp(b_0 + sex[t]*b_2 + b_1_years[years_i[t-1] + i]);
          // transition rate matrix, changes every time 
          Q[1,2] = q[1];
          Q[1,3] = 0;
          Q[1,4] = q[2];
          Q[2,1] = 0;
          Q[2,3] = q[3];
          Q[2,4] = q[4];
          Q[3,1] = 0;
          Q[3,2] = 0;
          Q[3,4] = q[5];
          Q[4,1] = 0;
          Q[4,2] = 0;
          Q[4,3] = 0;
          
          Q[1,1] = -sum(q[1:2]);
          Q[2,2] = -sum(q[3:4]);
          Q[3,3] = -q[5];
          Q[4,4] = 0;
          
          // first period over the previous age
          if(i == 1){ // tpm between previous year and next integer year
            gamma_temp = matrix_exp((years_i[t-1] + 1 - years[t-1])*Q);
          } else{
            if( i == per){ // tpm between last integer year and obs age
              gamma_temp = matrix_exp((years[t] - years_i[t])*Q);
            }else{ //full inte-observations integer years
              gamma_temp = matrix_exp(Q) ;
            }
          } // compute the tpm for all the inter-observations time
          gamma[t] = gamma[t] * gamma_temp ;
       }
    }
    
    // each row must sum to 1, compute log(ptm) between previous and current observation
    for(i in 1:N) {
      log_gamma[t][i] = log((gamma[t][i]+0.000000000001)/sum(0.000000000001+gamma[t][i]));
    }
    
    // transpose the log(tpms),each row corresponds to the transition into the state 
    log_gamma_tr[t] = log_gamma[t]';
    
    // account for the fact that dead is known without error
    if(y[t] == 4){
      Q_t = Q';
      for(n in 1:3){
        log_gamma_tr[t, 4, n] = log_sum_exp(to_vector(log_gamma[t,n]) + log(to_vector(Q_t[4])));
       }
     }

  }
 
 // likelihood computation 
  for (t in 1:T) {
    // initialise forward variable if first obs of track
    if(t == 1 || ID[t] != ID[t-1]){
      // compute nu 
      nu =  delta' * gamma[t]; 
      for(i in 1:3) { // conditional probabilities for admitted initial states
        delta_nu[i] = log((nu[i])/sum(nu[1:3]));
      }
      delta_nu[4] = 0;
    
      // forward algorithm implementation
      // first observation
      logp = delta_nu;
      // update log-like  
      for (n in 1:3) {
        logptemp[n] = logp[n] + log(B1[t]*asmat[n]) + categorical_lpmf(y[t] | to_vector(P[n])); 
      }
      logptemp[4] = logp[4] + categorical_lpmf(y[t] | to_vector(P[4]));
    } else{
        // update log-like  
         for (n in 1:3) {
          logptemp[n] = log_sum_exp(to_vector(log_gamma_tr[t,n]) + logp); 
          // F-V observations only when not dead
          if(y[t]!=4){logptemp[n] = logptemp[n] + log(B1[t]*asmat[n]) + categorical_lpmf(y[t] | to_vector(P[n]));} 
        }
        logptemp[4] = log_sum_exp(to_vector(log_gamma_tr[t,4]) + logp);
        logptemp[4] = logptemp[4] + categorical_lpmf(y[t] | to_vector(P[4]));
      }
    
    
    logp = logptemp;
    
    // add log forward variable to target at the end of each track
    if(t == T || ID[t+1] != ID[t]){
    target += log_sum_exp(logp);
    }
  }
}
