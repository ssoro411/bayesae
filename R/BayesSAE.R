library(rstan)
library(loo)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


basepath <- file.path(Sys.getenv("PROJECTS"), "bayesae")
stanpath <- file.path(basepath, "R/StanCodes")



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

    stanfit <- stan(file = file.path(stanpath, model_name), data = dat,
                    iter = iter, pars = c("sigma_sq", "beta", "theta",
                                          "log_lik"),
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

