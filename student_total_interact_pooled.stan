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
  real d_mu;
  real <lower = 0>d_sig;
  vector[8] q;
  real<lower = 0>sig;

}

model {
  // priors
 sig ~ exponential( 1 );
 q ~ normal( 0 , 1 );
 d_sig ~ exponential( 1 );
 d_mu ~ normal( 0 , 0.2 );
 c ~ normal( 0 , 0.2 );
 b ~ normal( 0 , 0.2 );
 a ~ normal( 0 , 0.2 );
 nu ~ gamma( 2 ,0.1 );

  //likelihood
  for (i in 1:n) {
    z[i] ~ student_t(nu,a[group[i]]+b[conchx[i]] + c[sex[i]] + d_mu + q[int_term[i]]*d_sig, sig);
  }


}

generated quantities{
    vector[n] log_lik;
     vector[n] mu;
     vector[8] d;
     d=d_mu + q*d_sig;
    for ( i in 1:n ) {
        mu[i] = a[group[i]]+b[conchx[i]] + c[sex[i]] + d_mu + q[int_term[i]]*d_sig;
    }
    for ( i in 1:n ) log_lik[i] = student_t_lpdf( z[i] | nu,mu[i] , sig );
}
