#' Univariate hierarchical Bayes approach to small area estimation.
#'
#' Hierarchical Bayes Approach to Small Area Estimation interfacing with Stan.
#' @name BayesSAE
#' @param formular formula
#' @param data Data frame with direct estimate and auxiliary variables.
#' @param Di \code{m} vector with sampling variance.
#' @param domain Vector with Domain names.
#' @param model There are three possible models. "FH" for Fay-Herriot model, "CAR" for conditional auto-regressive model and "SAR" for simultaneous auto-regressive model.
#' @param W Spatial matrix. If \code{model}="SAR", rowsum should be 1.
#' @param range Range of eigenvalues of W (only used if \code{model}="CAR").
#' @param logit.trans If true, it transforms simulated theta values to inv.logit(theta).
#' @param pars Parameters to be monitored.
#' @param iter Total iteration.
#' @param warmup Warm up. Default is "\code{floor}(\code{iter}/2)".
#' @param chains Number of chains. Default is 4.
#' @param control See the \code{rstan} package document.
#' @param open.progress Progress of chiain will be presented if it is \code{TRUE}.
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


#########################################################################
##  Function that conducts posterior simulation using Stan
#########################################################################

BayesSAE <- function(formula, data = NULL , Di = NULL, domain = NULL,
                        model = "FH", W = NULL, range = NULL, logit.trans=TRUE,
                        pars=c("sigma_sq", "beta", "theta","log_lik"),
                        iter = 1000, warmup = floor(iter/2), chains = 4,
                        control = list(max_treedepth=12, adapt_delta = 0.95),
                        open.progress = TRUE){

    this.call <- as.list( sys.call() )
    mf_ <- model.frame(formula, data = data)
    Y <- model.extract(mf_, "response")
    X <- model.matrix(formula, data = data)
    aux <- names(mf_)[-1]
    model_name <- paste("Stan_",model,".stan",sep="")

    if( is.null(domain ) ){ domain = rownames(data) }

    if(model == "FH"){
    dat <- list(m=dim(X)[1], p=dim(X)[2], y=Y, X=X, sDi= sqrt(Di) )
    }else if (model == "SAR"){
    dat <- list(m=dim(X)[1], p=dim(X)[2], y=Y, X=X, sDi= sqrt(Di), W=W)
    }else if (model == "CAR"){
    dat <- list(m=dim(X)[1], p=dim(X)[2], y=Y, X=X, sDi= sqrt(Di), W=W, rupper=range[2],rlower=range[1] )
     }


    stanfit <- stan(model_code = Model(model), model_name = model,
                    data = dat, pars = pars,
                    iter = iter, warmup = warmup, chains = chains,
                    open_progress = open.progress,
                    control = control )

    theta.smpl <- extract(stanfit, pars = "theta", permuted = FALSE)

    ll = extract_log_lik(stanfit)
    model_qual = list( LOO = loo(ll), WAIC = waic(ll) )

    if(logit.trans) theta.smpl <- expit(theta.smpl)
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




