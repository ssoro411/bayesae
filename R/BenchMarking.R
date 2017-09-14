#' Benchmarked Bayes Estimator
#'
#' Constrained Bayes estimator which minimizes posterior weighted squared error loss.
#' @param theta_b Bayes estimator.
#' @param w Weight vector associates with the benchmark.
#' @param t Pre-specified number.
#' @param phi Weight vector associates with weighted squared sum of square loss.
#' @param lambda Penalty parameter. When it is NULL, lambda = infinity and it will give exact benchmarking.
#' @return Resulting Benchmarked Bayes Estimator.
#' @export
#' @references
#' \insertRef{datta2011bench}{bayesae}


## Benchmarked Bayes Estimator

bbm = function(theta_b, w, t, phi, lambda=NA){
r = w/phi
s = sum(w*r)

if( is.na(lambda) ){
  theta_BM = theta_b + ( (s )^(-1) )*(t - sum(w*theta_b) )*r
  } else{
  theta_BM = theta_b + ( ( lambda/( s*lambda + 1 ) )^(-1) )*(t - sum(w*theta_b) )*r}

return(theta_BM)
}


