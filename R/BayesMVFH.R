#' Multivariate Fay-Herriot Model interfacing with Stan
#'
#' Bayesian approach to Multivariate Fay-Herriot models.
#' @name BayesMVFH
#' @param direct Direct estimates.
#' @param aux Auxiliary variables.
#' @param Di Pre-specified number.
#' @param domain Domain names.
#' @param pars Parameters to be monitored.
#' @param iter Total iteration.
#' @param warmup Warm up. Default is "iter/2".
#' @param chains Number of chains. Default is 4.
#' @param control See the "rstan" document.
#' @param open.progress Progress of chiain will be presented if it is TRUE.
#'
#' @return Simulated posterior sample from the Stan.
#' @export
#' @references
#'
#' \insertRef{carpenter2016stan}{bayesae}
#'
#' \insertRef{guo2016rstan}{bayesae}



library(rstan)
library(loo)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

BayesMVFH <- function(direct= NULL, aux = NULL , Di = NULL, domain = NULL,
                      pars = NA, iter = 1000, warmup = floor(iter/2), chains = 4,
                     control = list(max_treedepth=12, adapt_delta = 0.95),
                     open.progress = TRUE){

  this.call  <- as.list( sys.call() )
  model_name <- paste("Multivariate_FH.stan")

  dat = list(m=dim(aux)[1], s=dim(aux)[2],  p=dim(aux)[3], yi=direct, X=aux, Di= Di  )

  if( is.null(domain ) ){ domain = 1:dim(aux)[1] }

  stanfit <- stan(model_code = Model("MV"), model_name = "Multivariate_FH",
                  data = dat, pars = pars,
                  iter = iter, warmup = warmup, chains = chains,
                  open_progress = open.progress,
                  control = control )

  theta.smpl <- extract(stanfit, pars = "theta", permuted = FALSE)

  ll = extract_log_lik(stanfit)
  model_qual = list( LOO = loo(ll), WAIC = waic(ll) )

  posterior.summary <- data.frame(domain = domain,
                                  monitor(theta.smpl, digits_summary = 5,
                                          warmup = 0,
                                          probs = c(0.025, 0.50, 0.975),
                                          print = FALSE))[,-3]
  names(posterior.summary) <- c("domain", "mean", "sd",
                                "Q.025", "median", "Q.975", "n_eff",
                                "Rhat")
  stanfit.slots <- sapply(slotNames(stanfit), slot, object = stanfit,
                          simplify=F)

  result <- do.call("new", append(list("stanfit.sae",
                                       estimates = posterior.summary,
                                       fitness   = model_qual,
                                       model.call = this.call),
                                  stanfit.slots))
  return( result )
}




