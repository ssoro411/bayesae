#' Stan models
#'
#' Stan models for Fay-Herriot, simultaneous autoregressive and conditional autoregressive models.
#' @name Model_f
#' @param model Name of a model (\code{"FH"}, \code{"CAR"}, \code{"SAR"}, \code{"MVFH"}).
#' @import rstan loo boot
#' @export
#' @references
#'

#########################################################################
##  Function specifying selected model.
#########################################################################
Model_f = function(model) {
  if(model=="FH"){
    #########################################################################
    ##  Stan Model : Fay-Herriot
    #########################################################################
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
    model {
    for( i in 1:m)
    theta[i] ~ normal(X[i]'*beta, sqrt(sigma_sq));
    for( i in 1:m)
    y[i]     ~ normal(theta[i], sDi[i]);
    }
    generated quantities{
    vector[m] log_lik;            // This is for model cheking criteria (loo, WAIC)
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
      matrix[m,m]   I;           // Identity matrix
      }
      parameters {
      real<lower=-0.99999, upper=0.99999> rho;   // Spatial parameter
      real<lower=0>  sigma_sq;                   // Scale parameter
      vector[7]          beta;                   // Regression parameter
      vector[m]         theta;                   // Area characteristics
      }
      model {
      theta    ~ multi_normal_prec(X*beta, (1/sigma_sq)*((I-rho*(W'))*(I-rho*(W))) );
      for( i in 1:m)
      y[i]     ~ normal(theta[i], sDi[i]);
      }
      generated quantities{
      vector[m] log_lik;            // This is for model cheking criteria (loo, WAIC)
      for( i in 1:m)
      log_lik[i] = normal_lpdf(y[i]|theta[i],sDi[i]);
      }

      "
    } else if (model == "CAR") {
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
      model {
      theta    ~ multi_normal_prec(X*beta, (1/sigma_sq)*( (I-rho*(W)) ) );
      for( i in 1:m)
      y[i]     ~ normal(theta[i], sDi[i]);
      }
      generated quantities{
      vector[m] log_lik;            // This is for model cheking criteria (loo, WAIC)
      for( i in 1:m)
      log_lik[i] = normal_lpdf(y[i]|theta[i],sDi[i]);
      }

      "
    } else {
      m <-"
      data {
      int<lower=1>       m;      // number of small areas
      int<lower=1>       s;      // Dimension of data
      int<lower=1>       p;      // Number of auxiliary variables
      matrix[s,p]        X[m];   // Auxiliary variables
      vector[s]          yi[m];  // Response vector
      vector<lower=0>[s] Di[m];  // Sampling errors (Variance)
      }
      parameters {
      vector[p] beta;             //
      vector[s] theta[m];         // m, s x 1 area characteristics
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
      yi[i]     ~ multi_normal(theta[i], diag_matrix(Di[i]));
      }
      generated quantities{
      vector[m] log_lik;
      for( i in 1:m)
      log_lik[i] = multi_normal_cholesky_lpdf(yi[i]|theta[i],diag_matrix(Di[i]));
      }

      "
    }
  return(m)
  }






