//
data {
  int<lower = 0> m; // dimension of observed data (# of biomarkers)
  int<lower = 0> k; // number of latent factors
  int<lower = 0> n; // number of participants
  int<lower = 0> n_obs; // number of observations
  //int<lower = 1, upper = 2> group[n]; //group variable
  vector[m] y[n];
  
}

transformed data {
  int<lower = 1> p; // number of nonzero lower triangular factor loadings
  vector[m] zeros;

  zeros = rep_vector(0, m);
  p = k * (m - k) + k * (k - 1) / 2;
  
}

parameters {
  vector<lower = 0>[k] beta_diag;
  vector[p] beta_lower_tri;
  vector<lower = 0>[m] sigma_y; // residual error
  real<lower = 0> sigma_L; // hyperparameter for loading matrix elements
 // matrix[2,2] beta;
 // vector<lower = 0>[2] sigma_beta;
  //corr_matrix[k] omega;
}

transformed parameters {
  matrix[m, k] L;
  cov_matrix[m] Sigma;
  matrix[m, m] L_Sigma;
  //matrix[n,k] z; // factor scores
  
  
  { // temporary scope to define loading matrix L
    int idx = 0;
    for (i in 1:m){
      for (j in (i + 1):k){
        L[i, j] = 0; // constrain upper tri elements to zero
      }
    }
    for (j in 1:k){
      L[j, j] = beta_diag[j];
      for (i in (j + 1):m){
        idx = idx + 1;
        L[i, j] = beta_lower_tri[idx];
      }
    }
  }
  
  Sigma = multiply_lower_tri_self_transpose(L);
  for (i in 1:m) Sigma[i, i] = Sigma[i, i] + sigma_y[i];
  L_Sigma = cholesky_decompose(Sigma);
  
  //z scores
  //for (i in 1:n) {
//for (j in 1:k) {
//z[i, j] = L[, j]' * inverse(Sigma) * y[i,];
 // }
 // } 
}

model {
 
  //priors for model 1
  beta_lower_tri ~ normal(0, sigma_L);
  sigma_L ~ normal(0,0.4);
  sigma_y ~ normal(0,0.4);
  
  // priors for diagonal entries (Leung and Drton 2016)
  for (i in 1:k) {
    target += (k - i) * log(beta_diag[i]) - .5 * beta_diag[i] ^ 2 / sigma_L;
  }
  
  // likelihood for factor model
  //y ~ multi_normal_cholesky(zeros,L_Sigma);
  
}
