\name{BayesSAE}
\alias{BayesSAE}
\title{Hierarchical Bayes approach to small area estimation interfacing with Stan}
\description{
Posterior simulation conducted by Stan
}
\usage{
BayesSAE <- function(formula, data = NULL , Di = NULL, domain = NULL,
                        model = "FH", W = NULL, range = NULL, logit.trans=TRUE,
                        pars=c("sigma_sq", "beta", "theta","log_lik"),
                        iter = 1000, warmup = floor(iter/2), chains = 4,
                        control = list(max_treedepth=12, adapt_delta = 0.95),
                        open.progress = TRUE, ...)
}
\argumnets{
\item{formula}{formula}
\item{data}{Dataframe with direct estimate and auxiliary variables.}
\item{Di}{Sampling variance}
\item{domain}{}
\item{model}{"FH","CAR","SAR"}
\item{W}{Spatial matrix with 0 and 1.}
\item{range}{Range of eigenvalues of W (only used \item{model}="CAR")}
}}
}
\examples{
hello()
}
