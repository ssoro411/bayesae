\name{BayesSAE}
\alias{BayesSAE}
\title{Hierarchical Bayes approach to small area estimation interfacing with Stan}
\description{
Posterior simulation via Stan.
}
\usage{
BayesSAE <- function(formula, data = NULL , Di = NULL, domain = NULL,
                        model = "FH", W = NULL, range = NULL, logit.trans=TRUE,
                        pars=c("sigma_sq", "beta", "theta","log_lik"),
                        iter = 1000, warmup = floor(iter/2), chains = 4,
                        control = list(max_treedepth=12, adapt_delta = 0.95),
                        open.progress = TRUE, ...)
}
\arguments{
\item{formula}{formula}
\item{data}{Data frame with direct estimate and auxiliary variables.}
\item{Di}{Sampling variance.}
\item{domain}{Vector with Domain names.}
\item{model}{There are three possible models. "FH" for Fay-Herriot model, "CAR" for conditional auto-regressive model and "SAR" for simultaneous auto-regressive model.}
\item{W}{Spatial matrix. If \code{model}="SAR", rowsum should be 1.}
\item{range}{Range of eigenvalues of W (only used \code{model}="CAR")}
\item{logit.trans}{If true, it transforms simulated theta values to inv.logit(theta).}
\item{pars}{Parameters to be monitored.}
}
\examples{
To be updated.
}


