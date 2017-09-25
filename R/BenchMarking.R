#' Benchmarked bayes estimator
#'
#' Constrained Bayes estimator which minimizes posterior weighted squared error loss.
#' @name getBench
#' @param theta_b Bayes estimator.
#' @param t Pre-specified value to be benchmarked.
#' @param interval 95\% credible interval of \code{t}. If supplied, standard error of \code{t} is approximated by \code{range(interval)/4}.
#' @param w Weight vector associates with the benchmark.
#' @param phi Weight vector associates with weighted squared sum of square loss.
#' @param lambda Penalty parameter. Default is \eqn{\lambda = \infty} and it will conduct exact benchmark.
#' @return Resulting Benchmarked Bayes Estimator.
#' @import rstan loo boot
#' @export
#' @references
#'
#' \insertRef{datta2011bench}{bayesae}


## Benchmarked Bayes Estimator

getBench = function(theta_b, t, interval=NULL, w, phi, lambda = NULL){

r = w/phi
s = sum(w*r)

if( is.null(lambda) & is.null(interval) ){
  lambda = Inf
} else if (is.null(lambda) & is.null(interval)!=TRUE ){
  range  <- diff(interval)
  vari   <- as.vector( (range/4)^2 )
  lambda <- 1/vari
} else { lambda = lambda
         print("Both 95% credible interval and lambda are supplied. Only lambda is used")
}


if( is.finite(lambda) ){
  theta_BM = theta_b + ( lambda/(s*lambda+1) )*(t - sum(w*theta_b) )*r
  print("Lambda",lambda)
  } else{
  theta_BM = theta_b + ( (s )^(-1) )*(t - sum(w*theta_b) )*r
  print("Exact benchmark")
  }

return(list(Benchmarked_estimate = theta_BM, lambda= lambda ) )
}

