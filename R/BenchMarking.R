#' Benchmarked bayes estimator
#'
#' Constrained Bayes estimator which minimizes posterior weighted squared error loss.
#' @param theta_b Bayes estimator.
#' @param w Weight vector associates with the benchmark.
#' @param t Scalar or vector to be benchmarked.
#' @param phi Weight vector associates with weighted squared sum of square loss.
#' @param lambda Penalty parameter. Default is \eqn{\lambda = \infty} and it will conduct exact benchmark.
#' @return Resulting Benchmarked Bayes Estimator.
#' @import rstan loo
#' @export
#' @references
#'
#' \insertRef{datta2011bench}{bayesae}


## Benchmarked Bayes Estimator

bbm = function(theta_b, w, t, phi, lambda = Inf){
r = w/phi
s = sum(w*r)

if( is.finite(lambda) ){
  theta_BM = theta_b + ( lambda/(s*lambda+1) )*(t - sum(w*theta_b) )*r
  } else{
  theta_BM = theta_b + ( (s )^(-1) )*(t - sum(w*theta_b) )*r
  }

return(theta_BM)
}

