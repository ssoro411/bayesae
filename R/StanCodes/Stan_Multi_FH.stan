data {
  int<lower=1>       m;      // number of small areas
  int<lower=1>       s;      // Dimension of data
  int<lower=1>       p;      // Number of auxiliary variables
  matrix[s,p]        X[m];   // Auxiliary variables
  vector[s]          yi[m];  // Response vector
  vector<lower=0>[s] Si[m];  // Sampling errors (Known)
}
parameters { 
  vector[p] beta;             //
  vector[s] theta[m];         //
  cholesky_factor_cov[s]  L;  //
}                              
transformed parameters {
    vector[s] mu[m];
    for(i in 1:m)
    mu[i] = X[i]*beta;
}
model {
  for( i in 1:m)
  theta[i] ~ multi_normal_cholesky(mu[i], L);
  for( i in 1:m)
  yi[i]     ~ multi_normal(theta[i], diag_matrix(Si[i]));
}
generated quantities{
    vector[m] log_lik;
for( i in 1:m)
    log_lik[i] = multi_normal_cholesky_lpdf(yi[i]|theta[i],diag_matrix(Si[i]));
}
