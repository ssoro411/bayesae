#########################################################################
##  Function specifying selected model.
#########################################################################
Model = function(model) {
  if(model=="FH"){
    m <- "
    data {
    int<lower=0>       m;    // Number of small areas
    int<lower=0>       p;    // Number of auxiliary variables
    real            y[m];    // Direct Estimate
    real<lower=0> sDi[m];    // Square root of sampling variances
    vector[p]        X[m];   // Auxiliary variable matrix
    }
    parameters {
    real<lower=0>     sigma_sq;    // Scale parameter
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
    theta[i] ~ normal(mu[i], sqrt(sigma_sq));
    for( i in 1:m)
    y[i]     ~ normal(theta[i], sDi[i]);
    }
    generated quantities{
    vector[m] log_lik;
    for( i in 1:m)
    log_lik[i] = normal_lpdf(y[i]|theta[i],sDi[i]);
    }

    "} else if (model == "SAR") {
      #########################################################################
      ##  Stan Model : SAR
      #########################################################################

      m <- "
      data {
      int<lower=0>     m;        // Number of small areas
      int<lower=0>     p;        // Number of auxiliary variables
      real           y[m];       // Direct Estimate
      real<lower=0> sDi[m];     // Square root of sampling variances
      matrix<lower=0, upper=1>[m,m] W;           // Spatial matrix
      matrix[m,p]   X;           // Auxiliary variable matrix
      }
      parameters {
      real<lower=-0.99999, upper=0.99999> rho;   // Spatial parameter
      real<lower=0>  sigma_sq;                   // Scale parameter
      vector[7]          beta;                   // Regression parameter
      vector[m]         theta;                   // Area characteristics
      }
      transformed parameters {
      vector[m] mu;
      vector[m] ones;
      matrix<lower=0>[m,m] I;
      mu   = X*beta;
      ones = rep_vector(1,m);
      I    = diag_matrix(ones);
      }
      model {
      theta    ~ multi_normal_prec(mu, (1/sigma_sq)*((I-rho*(W))*(I-rho*(W'))) );
      for( i in 1:m)
      y[i]     ~ normal(theta[i], sDi[i]);
      }
      generated quantities{
      vector[m] log_lik;
      for( i in 1:m)
      log_lik[i] = normal_lpdf(y[i]|theta[i],sDi[i]);
      }

      "
    } else {
      #########################################################################
      ##  Stan Model : CAR
      #########################################################################

      m <- "
      data {
      int<lower=0>     m;        // Number of small areas
      int<lower=0>     p;        // Number of auxiliary variables
      real           y[m];       // Direct estimate
      real<lower=0> sDi[m];      // Square root of sampling variances
      matrix<lower=0, upper=1>[m,m]   W;           // Spatial matrix
      matrix[m,p]   X;           // Auxiliary variable matrix
      matrix[m,m]   I;           // Identity matrix
      real rupper;                   // Upper bound of rho
      real rlower;                   // Lower bound of rho
      }
      parameters {
      real<lower=rlower, upper=rupper> rho;  // Spatial parameter
      real<lower=0>  sigma_sq;               // Scale parameter
      vector[p]          beta;               // Regression parameter
      vector[m]         theta;               // Area characteristics
      }
      transformed parameters {
      vector[m] mu;
      mu = X*beta;
      }
      model {
      theta    ~ multi_normal_prec(mu, (1/sigma_sq)*( (I-rho*(W)) ) );
      for( i in 1:m)
      y[i]     ~ normal(theta[i], sDi[i]);
      }
      generated quantities{
      vector[m] log_lik;
      for( i in 1:m)
      log_lik[i] = normal_lpdf(y[i]|theta[i],sDi[i]);
      }

      "
}
  return(m)
  }






