#' Multivariate Fay-Herriot Model interfacing with Stan
#'
#' Bayesian approach to \code{s}-variate Fay-Herriot models using \code{stan}.
#' @name BayesMVFH
#' @param direct \code{m} by \code{s} matrix with direct estimate.
#' @param aux Array of auxiliary variable with \code{dim=c(m,s,p)}.
#' @param Di \code{m} by \code{s} matrix with sampling variance.
#' @param domain \code{m} vector with domain names.
#' @param pars Parameters to be monitored.
#' @param iter Total iteration.
#' @param warmup Warm up. Default is \code{floor}(\code{iter/2}).
#' @param chains Number of chains. Default is 4.
#' @param control See the \code{rstan} package document.
#' @param open.progress Progress of chiain will be presented if it is \code{TRUE}.
#' @import rstan loo
#' @return Simulated posterior sample from the \code{rstan}.
#' @export
#' @references
#'
#' \insertRef{carpenter2016stan}{bayesae}
#'
#' \insertRef{guo2016rstan}{bayesae}
#'
#' \insertRef{vehtari2014waic}{bayesae}

library(rstan)
library(loo)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


#########################################################################
setClass(
  Class = "stanfit.sae",
  slots = c(estimates = "data.frame", fitness = "list", model.call = "list"),
  contains = "stanfit"
)

#########################################################################
##  S4 method for extraction of SAEs from class stanfit.sae
#########################################################################

setGeneric("getEstimates",
           def = function(object){standardGeneric("getEstimates")})
setMethod("getEstimates", signature = "stanfit.sae",
          definition = function(object){object@estimates})



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

  posterior.summary <- data.frame(domain = domain, direct = as.vector(direct),
                                  monitor(theta.smpl, digits_summary = 5,
                                          warmup = 0,
                                          probs = c(0.025, 0.50, 0.975),
                                          print = FALSE))
  posterior.summary <- posterior.summary[,setdiff( colnames(posterior.summary), c("se_mean","n_eff","Rhat") )]
  posterior.summary <- posterior.summary[,c(1,2,3,4,6,5,7)]
  names(posterior.summary) <- c("domain","direct_est", "post_mean", "post_sd",
                                "median", "Q.025", "Q.975")
  stanfit.slots <- sapply(slotNames(stanfit), slot, object = stanfit,
                          simplify=F)

  result <- do.call("new", append(list("stanfit.sae",
                                       estimates = posterior.summary,
                                       fitness   = model_qual,
                                       model.call = this.call),
                                       stanfit.slots))
  return( result )
}




