data {
  int<lower=0>     m;        // Number of small areas
  int<lower=0>     p;        // Number of auxiliary variables
  real           y[m];       // Direct Estimate
  real<lower=0> sDi[m];      // Sampling Error
  matrix[m,m]   W;           // Spatial matrix
  matrix[m,p]   X;           // Auxiliary variable
  matrix[m,m]   I;           // Identity matrix
}
parameters {
  real<lower=-0.99999, upper=0.99999> rho;   // Spatial parameter 
  real<lower=0>  sigma_sq;                   // Scale parameter
  vector[7]          beta;                   // Regression parameter
  vector[m]         theta;                   // Area characteristics
}
transformed parameters {
  vector[m] mu;
  cov_matrix[m] G;
  mu = X*beta;
  G = (1/sigma_sq)*( (I-rho*(W))*(I-rho*(W'))  ) ;// SAR precision matrix
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
