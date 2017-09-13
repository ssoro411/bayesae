#' Benchmarked Bayes Estimator
#'
#' Constrained Bayes estimator which minimizes posterior weighted squared error loss.
#' @param theta_b Bayes estimator.
#' @param w Given weight vector.
#' @param t Pre-specified number.
#' @param phi
#' @param lambda
#' @return T square of the input
#' @export
#' @references Datta, Gauri Sankar, et al. "Bayesian benchmarking with applications to small area estimation." Test 20.3 (2011): 574-588.

## Benchmarked Bayes Estimator

BM = function(theta_b, w, t, phi, lambda=NA){
r = w/phi
s = sum(w*r)

if( is.na(lambda) ){
  theta_BM = theta_b + ( (s )^(-1) )*(t - sum(w*theta_b) )*r
  } else{
  theta_BM = theta_b + ( (s + lambda^(-1))^(-1) )*(t - sum(w*theta_b) )*r}

return(theta_BM)
}
