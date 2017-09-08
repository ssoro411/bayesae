data {
  int<lower=0>       m;    // Number of small areas
  int<lower=0>       p;    // Number of auxiliary variables
  real            y[m];    // Direct Estimate
  real<lower=0> sDi[m];    // Square root of sampling variances
  vector[p]        X[m];    // Auxiliary variable
}
parameters {
  real<lower=0>     sigma;       // Scale parameter
  vector[p]          beta;       // Regression parameter
  vector[m]         theta;       // Area characteristics
}
transformed parameters {
  vector[m] mu;
    for( i in 1:m)
    mu[i] = X[i]'*beta;
}
model {
  for( i in 1:m)
  theta[i] ~ normal(mu[i], sigma);
  for( i in 1:m)
  y[i]     ~ normal(theta[i], sDi[i]);
}
generated quantities{
vector[m] log_lik;
for( i in 1:m)
log_lik[i] = normal_lpdf(y[i]|theta[i],sDi[i]);
}
