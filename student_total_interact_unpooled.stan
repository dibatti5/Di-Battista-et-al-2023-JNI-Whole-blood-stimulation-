//
data {
  int<lower = 0> n; // number of participants
  int<lower = 1, upper = 2> group[n]; //group variable
  int<lower = 1, upper =2> conchx[n];
  int<lower = 1, upper =2> sex[n];
  int<lower=1,upper = 8> int_term[n];
  vector [n]z;
}

parameters {
  real<lower = 0> nu;
  vector [2]a;
  vector [2]b;
  vector [2]c;
  vector[8] d;
  real<lower = 0>sig;

}

model {
  // priors
 d ~ normal( 0 , 0.2 );
 c ~ normal( 0 , 0.2 );
 b ~ normal( 0 , 0.2 );
 a ~ normal( 0 , 0.2 );
 nu ~ gamma( 2 ,0.1 );
 sig ~ exponential( 1 );

  //likelihood
  for (i in 1:n) {
    z[i] ~ student_t(nu,a[group[i]]+b[conchx[i]] + c[sex[i]] + d[int_term[i]], sig);
  }


}

generated quantities{
    vector[n] log_lik;
     vector[n] mu;
    for ( i in 1:n ) {
        mu[i] = a[group[i]]+b[conchx[i]] + c[sex[i]] + d[int_term[i]];
    }
    for ( i in 1:n ) log_lik[i] = student_t_lpdf( z[i] | nu,mu[i] , sig );
}
