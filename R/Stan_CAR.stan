data {
  int<lower=0>     m;        // Number of small areas
  int<lower=0>     p;        // Number of auxiliary variables
  real           y[m];       // Direct estimate
  real<lower=0> sDi[m];      // Sampling standard deviation
  matrix[m,m]   W;           // Spatial matrix
  matrix[m,p]   X;           // Auxiliary variable
  matrix[m,m]   I;           // Identity matrix
  real rupper;                   // Upper bound of rho
  real rlower;                   // Lower bound of rho
}
parameters {
  real<lower=rlower, upper=rupper> rho;  // Spatial parameter
  real<lower=0>  sigma_sq;       // Scale parameter
  vector[p]          beta;       // Regression parameter
  vector[m]         theta;       // Characteristic of areas
}
transformed parameters {
  vector[m] mu;
  cov_matrix[m] G;
  mu = X*beta;
  G = (1/sigma_sq)*( (I-rho*(W)) ) ; // CAR precision matrix
}
model {
  theta    ~ multi_normal_prec(mu, G);
  for( i in 1:m)
  y[i]     ~ normal(theta[i], sDi[i]);
}
generated quantities{
    vector[m] log_lik;
    for( i in 1:m)
    log_lik[i] = normal_lpdf(y[i]|theta[i],sDi[i]);
}

