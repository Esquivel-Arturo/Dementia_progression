data {
  int<lower=0> N; // number of states
  int<lower=0> R; // number of rates modelled 
  int<lower=1> T; // length of data set  
  int ID[T]; // id's 
  real pib[T]; // observations (log(PIB-1))
  real thick[T]; // observations
  real mmse[T]; // observations
  int dem[T]; // observations
  int sex[T]; // Covariates
  real age[T]; // age
  int age_i[T]; // Integer age
  real educ[T]; // years of education 
  int apoe[T]; // apoe4 precense
  real ntest[T]; // number of mmse taken
  int iyears[T]; // number of years enrolled
  matrix[71, 8] B; // splines
  int obstype[T]; // indicator of death observed
  int missp[T]; // indicator that pib is not missing
  int misst[T]; // indicator that thick is not missing
  
}  

parameters {
  vector<lower=0, upper=1>[2] p; // misclassification probabilitites
  simplex[4] del; // initial distribution 
  real b_0nD_Dead; // Intercepts
  real b_02 ;
  real<lower=b_02> b_04;
  real b_06;
  real b_07 ;
  real<lower=b_07> b_09;
  real<lower=b_06> b_011;
  real<lower=b_0nD_Dead> b_012;
  real<lower=b_012> b_013;
  real<lower=0> b_1nD_Dead; // coefficients, nD_Dead are the non-dementia to death coefficients
  real<lower=0> b_12 ;
  real<lower=0> b_14;
  real<lower=0> b_16;
  real<lower=0> b_17 ;
  real<lower=0> b_19;
  real<lower=0> b_111;
  real<lower=0> b_112;
  real<lower=0> b_113;
  real b_21;
  real b_2nD_Dead;
  real b_22 ;
  real b_24;
  real b_26;
  real b_27 ;
  real b_29;
  real b_211;
  real b_212;
  real b_213;
  real b_31;
  real b_3nD_Dead;
  real b_32 ;
  real b_34;
  real b_36;
  real b_37 ;
  real b_39;
  real b_311;
  real b_312;
  real b_313;
  real b_41;
  real b_4nD_Dead;
  real b_42 ;
  real b_44;
  real b_46;
  real b_47 ;
  real b_49;
  real b_411;
  real b_412;
  real b_413;
  real<upper=0> c; // death rate bias
  real<lower=0> d; // adjustment per year enrolled
  vector[8] sp; // spline coefs
  ordered[2] mu_a; // level dependent means
  ordered[2] mu_n;
  real<lower=0> sigma_p; // pib sd
  real<lower=0> sigma_t; // thicknes sd
  vector[5] alpha711; // mmse mean coefficients 
  real alpha1;
  vector<upper=alpha1>[2] alpha23;
  real<upper=min(alpha23)> alpha4;
  real<upper=alpha4> alpha5;
  real<upper=alpha5> alpha6;
  real<lower=0> sigma_m; // mmse sd
}

model {
  vector[N] logp; // log marginal likelihood
  matrix[N, N] Q; //Transition rates matrix, changes for each obs
  matrix[N, N] Q_t; //Auxiliary matrix for transpose of Q
  vector[N] logptemp; // auxiliary variable to compute the likelihood
  matrix[N,N] gamma[T]; // Transition probability matrix between t-th obs and previous obs
  matrix[N,N] log_gamma[T]; // log of the ptm
  matrix[N,N] log_gamma_tr[T]; // transposed 
  vector[R] q; // transition rates 
  int per; // auxiliary variable for the number of integer years in each transition
  matrix[N,N] gamma_temp; // auxiliary matrix (computed each integer year) to build ptm 
  vector[N] I; // auxiliary vector to build identity matrix
  row_vector[N] nu; // probabilities for first observation
  vector[N] delta_nu; // log probabilities of modified initial distribution
  vector[N] delta; //auxiliary for initial distribution 
  vector[2] g; // vector for g
  vector[6] mu; //vector for the means' intercepts of the mmse emission distribution
  vector[11] alpha; //vector with all mmse coefficients

  // priors for betas
  b_0nD_Dead ~ normal(-4.41, 0.1); // non-dem to dead
  
  b_012 ~ normal(-4, 1); // dem to dead
  b_013 ~ normal(-4, 1);
  
   // rest
  b_02 ~ normal(-3, 1);
  b_04 ~ normal(-3, 1);
  b_06 ~ normal(-3, 1);
  b_07 ~ normal(-3, 1);
  b_09 ~ normal(-3, 1);
  b_011 ~ normal(-3, 1);
  
  b_1nD_Dead ~ normal(0.094, 0.01); // non-dem to dead
  
  b_112 ~ normal(0.1, 0.05); // dem to dead
  b_113 ~ normal(0.1, 0.05);
  
  b_12 ~ normal(0.1, 0.05); // rest
  b_14 ~ normal(0.1, 0.05);
  b_16 ~ normal(0.1, 0.05);
  b_17 ~ normal(0.1, 0.05);
  b_19 ~ normal(0.1, 0.05);
  b_111 ~ normal(0.1, 0.05);
  
  b_2nD_Dead  ~ normal(0.47, 0.05); // non-dem to dead
  
  b_212 ~ normal(0, 1); // dem to dead
  b_213 ~ normal(0, 1);
  
  b_21 ~ normal(0, 1); // rest
  b_22 ~ normal(0, 1);
  b_24 ~ normal(0, 1);
  b_26 ~ normal(0, 1);
  b_27 ~ normal(0, 1);
  b_29 ~ normal(0, 1);
  b_211 ~ normal(0, 1);
  
  b_3nD_Dead  ~ normal(0, 0.1);
  b_31  ~ normal(0, 0.1);
  b_32  ~ normal(0, 0.1);
  b_34  ~ normal(0, 0.1);
  b_36  ~ normal(0, 0.1);
  b_37  ~ normal(0, 0.1);
  b_39  ~ normal(0, 0.1);
  b_311  ~ normal(0, 0.1);
  b_312  ~ normal(0, 0.1);
  b_313  ~ normal(0, 0.1);
  
  b_4nD_Dead ~ normal(0, 1);
  b_41 ~ normal(0,1);
  b_42  ~ normal(0,1);
  b_44  ~ normal(0,1);
  b_46  ~ normal(0,1);
  b_47  ~ normal(0,1);
  b_49  ~ normal(0,1);
  b_411  ~ normal(0,1);
  b_412  ~ normal(0,1);
  b_413  ~ normal(0,1);
  
  // death rate bias priors
  c ~ normal(-0.75, 0.375);
  d ~ normal(0.15, 0.3);
  
  // spline coeffs priors
  sp[1] ~ normal(-5, 1);
  sp[2] ~ normal(-4, 2);
  sp[3] ~ normal(-3, 2);
  sp[4] ~ normal(-2, 2);
  sp[5] ~ normal(-1, 2);
  sp[6] ~ normal(0, 3);
  sp[7] ~ normal(1, 3);
  sp[8] ~ normal(2, 3);
  
  // PIB priors
  mu_a[1] ~ normal(-1.3, 0.2);
  mu_a[2] ~ normal(-0.5, 0.2);
  sigma_p ~ normal(0.1, 2);
  
  // Thickness priors
  mu_n[1] ~ normal(2.34, 0.2);
  mu_n[2] ~ normal(3.14, 0.2); // this one corresponds to low burden
  sigma_t ~ normal(0.1, 2);
  
  // MMSE priors
  alpha1 ~ normal(-0.28, 0.75);
  alpha23 ~ normal(-0.28, 0.75);
  alpha4 ~ normal(-0.28, 0.75);
  alpha5 ~ normal(-7.3, 3);
  alpha6 ~ normal(-7.3, 3);
  alpha711 ~ normal(0, 1);
  sigma_m ~ normal(0.5, 2);
  
  // misclassification priors
  p ~ normal(0.05, 0.25) ;
  
  // initial distribution 
  del[1] ~ normal(0.966, 0.1);
  del[2] ~ normal(0.03, 0.1);
  del[3] ~ normal(0.003, 0.1);
  del[4] ~ normal(0.003, 0.1);
  
  // assign all alphas to the same vector
  alpha[1] = alpha1;
  alpha[2] = alpha23[1];
  alpha[3] = alpha23[2];
  alpha[4] = alpha4;
  alpha[5] = alpha5;
  alpha[6] = alpha6;
  
  for(i in 1:5){
    alpha[6+i] = alpha711[i];
  }
  
  // give values to I
  for(i in 1:N){
    I[i] = 1;
  }
  // only the first 4 states are possible at baseline
  for(i in 1:4){
    delta[i] = del[i];
  }
  
  for(i in 5:7){
    delta[i] = 0;
  }
  
  // derive array of (log-)transition probabilities
  for(t in 1:T){
    // define g
    g[1] = c + d*iyears[t];
    g[2] = 0;
    // build NxN identity matrix
    gamma[t] = diag_matrix(I);
    // for first observation
    if( t == 1 || ID[t] != ID[t-1]) {
      // number of integer periods to consider between first obs time and baseline (50)
      per = age_i[t] + 19;
      for(i in 1:per){
        // transitions rates at each year in p 
        q[1] = exp(B[i] * sp + b_21*sex[t] + b_31*educ[t] + b_41*apoe[t]) ;
        q[2] = exp(b_02 + b_12 * (i - 19) + b_22*sex[t] + b_32*educ[t] + b_42*apoe[t]) ;
        q[3] = exp(b_0nD_Dead + b_1nD_Dead * (i - 19) + b_2nD_Dead*sex[t] + b_3nD_Dead*educ[t] + b_4nD_Dead*apoe[t] + min(g));
        q[4] = exp(b_04 + b_14 * (i - 19) + b_24*sex[t] + b_34*educ[t] + b_44*apoe[t]) ;
        q[5] = exp(b_0nD_Dead + b_1nD_Dead * (i - 19) + b_2nD_Dead*sex[t] + b_3nD_Dead*educ[t] + b_4nD_Dead*apoe[t] + min(g));
        q[6] = exp(b_06 + b_16 * (i - 19) + b_26*sex[t] + b_36*educ[t] + b_46*apoe[t]);
        q[7] = exp(b_07 + b_17 * (i - 19) + b_27*sex[t] + b_37*educ[t] + b_47*apoe[t]);
        q[8] = exp(b_0nD_Dead + b_1nD_Dead * (i - 19) + b_2nD_Dead*sex[t] + b_3nD_Dead*educ[t] + b_4nD_Dead*apoe[t] + min(g));
        q[9] = exp(b_09 + b_19 * (i - 19) + b_29*sex[t] + b_39*educ[t] + b_49*apoe[t]) ;
        q[10] = exp(b_0nD_Dead + b_1nD_Dead * (i - 19) + b_2nD_Dead*sex[t] + b_3nD_Dead*educ[t] + b_4nD_Dead*apoe[t] + min(g));
        q[11] = exp(b_011 + b_111 * (i - 19) + b_211*sex[t] + b_311*educ[t] + b_411*apoe[t]);
        q[12] = exp(b_012 + b_112 * (i - 19) + b_212*sex[t] + b_312*educ[t] + b_412*apoe[t]);
        q[13] = exp(b_013 + b_113 * (i - 19) + b_213*sex[t] + b_313*educ[t] + b_413*apoe[t]);
        
        // transition rate matrix, changes every time 
        Q[1,2] = q[1];
        Q[1,3] = q[2];
        Q[1,4] = 0;
        Q[1,5] = 0;
        Q[1,6] = 0;
        Q[1,7] = q[3];
        
        Q[2,1] = 0;
        Q[2,3] = 0;
        Q[2,4] = q[4];
        Q[2,5] = 0;
        Q[2,6] = 0;
        Q[2,7] = q[5];
        
        Q[3,1] = 0;
        Q[3,2] = 0;
        Q[3,4] = q[6];
        Q[3,5] = q[7];
        Q[3,6] = 0;
        Q[3,7] = q[8];
        
        Q[4,1] = 0;
        Q[4,2] = 0;
        Q[4,3] = 0;
        Q[4,5] = 0;
        Q[4,6] = q[9];
        Q[4,7] = q[10];
        
        Q[5,1] = 0;
        Q[5,2] = 0;
        Q[5,3] = 0;
        Q[5,4] = 0;
        Q[5,6] = q[11];
        Q[5,7] = q[12];
        
        Q[6,1] = 0;
        Q[6,2] = 0;
        Q[6,3] = 0;
        Q[6,4] = 0;
        Q[6,5] = 0;
        Q[6,7] = q[13];
        
        Q[7,1] = 0;
        Q[7,2] = 0;
        Q[7,3] = 0;
        Q[7,4] = 0;
        Q[7,5] = 0;
        Q[7,6] = 0;
        
  
        Q[1,1] = -sum(q[1:3]);
        Q[2,2] = -sum(q[4:5]);
        Q[3,3] = -sum(q[6:8]);
        Q[4,4] = -sum(q[9:10]);
        Q[5,5] = -sum(q[11:12]);
        Q[6,6] = -q[13];
        Q[7,7] = 0;
        
        // first period 
        if(i == 1){
          if(age[t]< -17){ // first observation within the first integer year
            gamma_temp = matrix_exp((age[t] + 18)*Q);
          } else{ // integer years within the inter-observations period
            gamma_temp = matrix_exp(Q);
          }
        } else{ // last inger year of the period
          if( i == per){
            gamma_temp = matrix_exp((age[t] - age_i[t])*Q);
          } else{
            gamma_temp = matrix_exp(Q) ;
          }
        } // compute the tpm for all the inter-observations time
        gamma[t] = gamma[t] * gamma_temp ;
      }
    } else{ // when it's not the first observation of the individual
        // number of integer periods to consider between current and previous observation
        per = age_i[t] - age_i[t-1] + 1;
        for(i in 1:per){ // transitions rates at each year in p 
           // transitions rates at each year in p 
          q[1] = exp(B[age_i[t] + 18 + i] * sp + b_21*sex[t] + b_31*educ[t] + b_41*apoe[t]) ;
          q[2] = exp(b_02 + b_12 * (age_i[t-1]+i-1) + b_22*sex[t] + b_32*educ[t] + b_42*apoe[t]) ;
          q[3] = exp(b_0nD_Dead + b_1nD_Dead * (age_i[t-1]+i-1) + b_2nD_Dead*sex[t] + b_3nD_Dead*educ[t] + b_4nD_Dead*apoe[t] + min(g));
          q[4] = exp(b_04 + b_14 * (age_i[t-1]+i-1) + b_24*sex[t] + b_34*educ[t] + b_44*apoe[t]) ;
          q[5] = exp(b_0nD_Dead + b_1nD_Dead * (age_i[t-1]+i-1) + b_2nD_Dead*sex[t] + b_3nD_Dead*educ[t] + b_4nD_Dead*apoe[t] + min(g));
          q[6] = exp(b_06 + b_16 * (age_i[t-1]+i-1) + b_26*sex[t] + b_36*educ[t] + b_46*apoe[t]);
          q[7] = exp(b_07 + b_17 * (age_i[t-1]+i-1) + b_27*sex[t] + b_37*educ[t] + b_47*apoe[t]);
          q[8] = exp(b_0nD_Dead + b_1nD_Dead * (age_i[t-1]+i-1) + b_2nD_Dead*sex[t] + b_3nD_Dead*educ[t] + b_4nD_Dead*apoe[t] + min(g));
          q[9] = exp(b_09 + b_19 * (age_i[t-1]+i-1) + b_29*sex[t] + b_39*educ[t] + b_49*apoe[t]) ;
          q[10] = exp(b_0nD_Dead + b_1nD_Dead * (age_i[t-1]+i-1) + b_2nD_Dead*sex[t] + b_3nD_Dead*educ[t] + b_4nD_Dead*apoe[t] + min(g));
          q[11] = exp(b_011 + b_111 * (age_i[t-1]+i-1) + b_211*sex[t] + b_311*educ[t] + b_411*apoe[t]);
          q[12] = exp(b_012 + b_112 * (age_i[t-1]+i-1) + b_212*sex[t] + b_312*educ[t] + b_412*apoe[t]);
          q[13] = exp(b_013 + b_113 * (age_i[t-1]+i-1) + b_213*sex[t] + b_313*educ[t] + b_413*apoe[t]);       
            // transition rate matrix, changes every time 
          Q[1,2] = q[1];
          Q[1,3] = q[2];
          Q[1,4] = 0;
          Q[1,5] = 0;
          Q[1,6] = 0;
          Q[1,7] = q[3];
        
          Q[2,1] = 0;
          Q[2,3] = 0;
          Q[2,4] = q[4];
          Q[2,5] = 0;
          Q[2,6] = 0;
          Q[2,7] = q[5];
        
          Q[3,1] = 0;
          Q[3,2] = 0;
          Q[3,4] = q[6];
          Q[3,5] = q[7];
          Q[3,6] = 0;
          Q[3,7] = q[8];
        
          Q[4,1] = 0;
          Q[4,2] = 0;
          Q[4,3] = 0;
          Q[4,5] = 0;
          Q[4,6] = q[9];
          Q[4,7] = q[10];
        
          Q[5,1] = 0;
          Q[5,2] = 0;
          Q[5,3] = 0;
          Q[5,4] = 0;
          Q[5,6] = q[11];
          Q[5,7] = q[12];
        
          Q[6,1] = 0;
          Q[6,2] = 0;
          Q[6,3] = 0;
          Q[6,4] = 0;
          Q[6,5] = 0;
          Q[6,7] = q[13];
        
          Q[7,1] = 0;
          Q[7,2] = 0;
          Q[7,3] = 0;
          Q[7,4] = 0;
          Q[7,5] = 0;
          Q[7,6] = 0;
        
  
          Q[1,1] = -sum(q[1:3]);
          Q[2,2] = -sum(q[4:5]);
          Q[3,3] = -sum(q[6:8]);
          Q[4,4] = -sum(q[9:10]);
          Q[5,5] = -sum(q[11:12]);
          Q[6,6] = -q[13];
          Q[7,7] = 0;
          
          // first period over the previous age
          if(i == 1){ // tpm between previous year and next integer year
            gamma_temp = matrix_exp((age_i[t-1] + 1 - age[t-1])*Q);
          } else{
            if( i == per){ // tpm between last integer year and obs age
              gamma_temp = matrix_exp((age[t] - age_i[t])*Q);
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
    
    // transpose the log(tpms), each row corresponds to the transition into the state 
    log_gamma_tr[t] = log_gamma[t]';
    
    // account for the fact that dead is known without error
     if(obstype[t] == 2){
       Q_t = Q';
       for(n in 1:6){
         log_gamma_tr[t, 7, n] = log_sum_exp(to_vector(log_gamma[t,n]) + log(to_vector(Q_t[7])));
        }
      }
}
  
 
 // likelihood computation 
  for(t in 1:T){
      // define mu 
      mu = alpha[1:6] + alpha[7]*age[t] + alpha[8]*sex[t] + alpha[9]*educ[t] + alpha[10]*apoe[t] + alpha[11]*ntest[t] ;

    // initialise forward variable if first obs of track
    if(t == 1 || ID[t] != ID[t-1]){
      // compute nu 
      nu =  delta' * gamma[t]; 
      for(i in 1:4) { // conditional probabilities for admitted initial states
        delta_nu[i] = log((nu[i]+0.0000000000001)/sum(nu[1:4]+0.0000000000001));
      }
      delta_nu[5] = 0;
      delta_nu[6] = 0;
      delta_nu[7] = 0;
    
      // forward algorithm implementation
      // first observation
      logp = delta_nu;
      
      // update log-like  
      for (n in 1:6) {
        if(n == 1 || n == 3 || n == 5){ // states that correspond to one PIB mean 
          logptemp[n] = logp[n] + normal_lpdf(pib[t] | mu_a[1], sigma_p)*missp[t]; 
        } else {
           // states that correspond to the other PIB mean 
          logptemp[n] = logp[n] + normal_lpdf(pib[t] | mu_a[2], sigma_p)*missp[t]; 
          }
        
        
        if(n == 1 || n == 2){ // states that correspond to one thickness mean 
          logptemp[n] = logptemp[n] + normal_lpdf(thick[t] | mu_n[2], sigma_t)*misst[t]; 
        } else {
           // states that correspond to the other thickness mean 
          logptemp[n] = logptemp[n] + normal_lpdf(thick[t] | mu_n[1], sigma_t)*misst[t]; 
          }          
        
         logptemp[n] = logptemp[n] + normal_lpdf(mmse[t] | mu[n], sigma_m);  
         
        if(n == 1 || n == 2 || n == 3 || n == 4) { // states that correspond to no dementia 
          logptemp[n] = logptemp[n] + bernoulli_lpmf(dem[t] | p[1]); 
        } else {
           // states that correspond to dementia
          logptemp[n] = logptemp[n] + bernoulli_lpmf(1 - dem[t] | p[2]); 
          }
        }
        
      logptemp[7] = logp[7];
      
    }else{
      // update log-like  
      for (n in 1:6) {
        logptemp[n] = log_sum_exp(to_vector(log_gamma_tr[t,n]) + logp);
        if(obstype[t] != 2){ 
          if(n == 1 || n == 3 || n == 5){ // states that correspond to one PIB mean 
            logptemp[n] = logptemp[n] + normal_lpdf(pib[t] | mu_a[1], sigma_p)*missp[t]; 
        } else {
          // states that correspond to the other PIB mean 
          logptemp[n] = logptemp[n] + normal_lpdf(pib[t] | mu_a[2], sigma_p)*missp[t]; 
          }
      
        
        if(n == 1 || n == 2){ // states that correspond to one thickness mean 
          logptemp[n] = logptemp[n] + normal_lpdf(thick[t] | mu_n[2], sigma_t)*misst[t]; 
        } else {
          // states that correspond to the other thickness mean 
          logptemp[n] = logptemp[n] + normal_lpdf(thick[t] | mu_n[1], sigma_t)*misst[t]; 
          }          
        
         logptemp[n] = logptemp[n] + normal_lpdf(mmse[t] | mu[n], sigma_m);  
        
        
        if(n == 1 || n == 2 || n == 3 || n == 4) { // states that correspond to one no dementia 
          logptemp[n] = logptemp[n] + bernoulli_lpmf(dem[t] | p[1]); 
        } else {
          // states that correspond to dementia
          logptemp[n] = logptemp[n] + bernoulli_lpmf(1 - dem[t] | p[2]); 
        }
      }
    }
      logptemp[7] = log_sum_exp(to_vector(log_gamma_tr[t,7]) + logp);
  }
    
    logp = logptemp;
    
    // add log forward variable to target at the end of each track
    if(t == T || ID[t+1] != ID[t]){
    target += log_sum_exp(logp);
    }
  }
}
