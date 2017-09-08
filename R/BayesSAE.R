library(rstan)
library(loo)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())





#########################################################################
setClass(
    Class = "stanfit.sae",
    slots = c(estimates = "data.frame", model.call = "list"),
    contains = "stanfit"
)
#########################################################################
##  S4 method for extraction of SAEs from class stanfit.sae
#########################################################################
setGeneric("getEstimates",
           def = function(object){standardGeneric("getEstimates")})
setMethod("getEstimates", signature = "stanfit.sae",
          definition = function(object){object@estimates})


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

"} else if (model == "SAR") {
#########################################################################
##  Stan Model : SAR
#########################################################################

m <- "
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

"
} else {
#########################################################################
##  Stan Model : SAR
#########################################################################

m <- "
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

"
}
return(m)
}







BayesSAE <- function(formula, data = NULL , Di = NULL, domain = NULL,
                        model = "FH", W = NULL, range = NULL,
                        iter = 2000, warmup = floor(iter/2), chains = 4,
                        open.progress = TRUE, ...){

    this.call <- as.list( sys.call() )
    mf_ <- model.frame(formula, data = data)
    Y <- model.extract(mf_, "response")
    X <- model.matrix(formula, data = data)
    aux <- names(mf_)[-1]
    model_name <- paste("Stan_",model,".stan",sep="")

    if(model == "FH"){
    dat <- list(m=dim(X)[1], p=dim(X)[2], y=Y, X=X, sDi= sqrt(Di) )
    }else if (model == "SAR"){
    dat <- list(m=dim(X)[1], p=dim(X)[2], y=Y, X=X, sDi= sqrt(Di), W=WSAR, I=I )
    }else {
    dat <- list(m=dim(X)[1], p=dim(X)[2], y=Y, X=X, sDi= sqrt(Di), W=WSAR, I=I, rupper=range[2],rlower=range[1] )
     }

    # stanfit <- stan(file = file.path(stanpath, model_name), data = dat,
    #                 iter = iter, pars = c("sigma_sq", "beta", "theta",
    #                                       "log_lik"),
    #                 warmup = warmup, chains = chains,
    #                open_progress=open.progress,
    #                control = list(max_treedepth=15, adapt_delta = 0.99))

    stanfit <- stan(model_code = Model(model), ,model_name = model,
                    data = dat, iter = iter, 
                    pars = c("sigma_sq", "beta", "theta","log_lik"),
                    warmup = warmup, chains = chains,
                    open_progress=open.progress,
                    control = list(max_treedepth=15, adapt_delta = 0.99))
    
    theta.smpl <- extract(stanfit, pars = "theta", permuted = FALSE)
    if(logit.trans) theta.smpl <- expit(theta.smpl)
    posterior.summary <- data.frame(domain = rownames(data),
                            monitor(theta.smpl, digits_summary = 5,
                                    probs = c(0.025, 0.50, 0.975),
                                    print = FALSE))
    names(posterior.summary) <- c("domain", "mean", "se_mean", "sd",
                                    "Q.025", "median", "Q.975", "n_eff",
                                  "Rhat")
    stanfit.slots <- sapply(slotNames(stanfit), slot, object = stanfit,
                            simplify=F)
    result <- do.call("new", append(list("stanfit.sae",
                                         estimates = posterior.summary,
                                         model.call = this.call),
                                    stanfit.slots))
    result
}


