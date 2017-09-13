#' Benchmarked Bayes Estimator
#'
#' Constrained Bayes estimator which minimizes posterior weighted squared error loss.
#' @param theta_b Bayes estimator.
#' @param w Given weight vector.
#' @param t Pre-specified number.
#' @param phi
#' @param lambda
#' @return Resulting Benchmarked Bayes Estimator.
#' @export
#' @references
#' \insertRef{Rpack:bibtex}{Rdpack}
#' \insertRef{datta2011bench}{bibtex}


## Benchmarked Bayes Estimator

bbm = function(theta_b, w, t, phi, lambda=NA){
r = w/phi
s = sum(w*r)

if( is.na(lambda) ){
  theta_BM = theta_b + ( (s )^(-1) )*(t - sum(w*theta_b) )*r
  } else{
  theta_BM = theta_b + ( (s + lambda^(-1))^(-1) )*(t - sum(w*theta_b) )*r}

return(theta_BM)
}


