//
data {
  int<lower = 0> m; // dimensions - number of biomarkers
  int<lower = 0> k; // number of latent factors
  int<lower = 0> n; // number of participants
  int<lower = 0> n_obs; // number of observations
  int<lower = 1, upper = n> row_obs[n_obs]; // row indices for observations
  int<lower = 1, upper = m> col_obs[n_obs]; // column indices for observations
  vector[n_obs] y_obs;
  int<lower=0> n_miss; // number of missing observations
  int<lower = 1, upper = n> row_miss[n_miss];
  int<lower = 1, upper = m> col_miss[n_miss];
  vector[n_miss] Lo;
  vector[n_miss] Up;
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
  vector<lower = 0, upper = 1>[n_miss] y_miss;
}

transformed parameters {
  matrix[m, k] L;
  cov_matrix[m] Sigma;
  matrix[m, m] L_Sigma;
  vector[m] y[n];
  vector[n_miss] y_miss2 = Lo + (Up - Lo) .* y_miss; //constrain imputed parameters between 0 and min sample value.
  matrix[n,k] z; // factor scores
  
  // build n by m matrix y by combining observed and missing observations
  for (i in 1:n_obs) {
    y[row_obs[i]][col_obs[i]] = y_obs[i];
   
  }
  for (i in 1:n_miss) {
    y[row_miss[i]][col_miss[i]] = y_miss2[i]; 
   
  }

  
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
  for (i in 1:n) {
for (j in 1:k) {
z[i, j] = L[, j]' * inverse(Sigma) * y[i,];
  }
  } 
}

model {
  // priors
  beta_lower_tri ~ normal(0, sigma_L);
  sigma_L ~ normal(0,0.5);
  sigma_y ~ normal(0,0.5);
  // priors for diagonal entries (Leung and Drton 2016)
  for (i in 1:k) {
    target += (k - i) * log(beta_diag[i]) - .5 * beta_diag[i] ^ 2 / sigma_L;
  }

  // likelihood
  y ~ multi_normal_cholesky(zeros, L_Sigma); 
}



