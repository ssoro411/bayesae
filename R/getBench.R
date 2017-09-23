#' Returning penalty parameter for benchmarked Bayes estimator based on supplied sample.
#'
#' Provide penalty parameter for Benchmarked Bayes Estimator.
#' @name getBench
#' @param fit \code{stan} ofr \code{array} of sample to be benchmarked. Target \eqn{t} will be the sample weighted average.
#' @param weight Weight vector associates with the benchmark.
#' @param par Parameter name to be benchmarked. Default is \eqn{\theta}.
#' @param interval 95\% credible interval of \eqn{t}. If supplied, other inputs are suppressed and standard error of \eqn{t} is approximated by \code{range(interval)/4}.
#' @export
#' @references
#'
#' \insertRef{datta2011bench}{bayesae}

# setwd( "/Users/heecheolchung/Dropbox/Research/SAE/SAE_Kenyea/MV_FH/MVFH_IndMH-12-30-16" )
# load("mv.RData")

# setwd( "/Users/heecheolchung/Dropbox/Research/SAE/saepackage" )
# load("Stan_Result3small.RData")

#########################################################################
# setClass(
    # Class = "bench.sae",
    # slots = c(result="list"),
            # )

#########################################################################
##  Function that gives target and penalty based on supplied sample.
#########################################################################

getBench <- function(fit, weight = NULL, par = "theta", interval = NULL ){

if( is.null(interval) ){
  if( class(fit)[1] == "stanfit" | class(fit)[1] == "stanfit.sae" ){
    postsam <- extract(fit, par="theta", permuted = TRUE)
  } else if ( class(fit)[1] == "array" | class(fit)[1] == "matrix"){
    postsam <- fit
  } else { stop(paste( "Data class should be an stanfit or array or matrix")) }

#str(postsam)
dims  <- dim(postsam$theta)
ndims <- length(dims)
benchsam <- postsam$theta
#str(benchsam)

# Supplied weight
if( is.null(weight) ){ weight <- rep(1,dims[2]) }

  tsample <- array( apply( benchsam, setdiff( 1:ndims,2) , function(x) sum(x*weight) ),
				            dim = c(dims[1],ifelse( ndims==3, dims[3], 1 )) )


  t       <- apply( tsample, 2, mean)
  vari    <- apply( tsample, 2, var )
  lambda  <- 1/vari
  intval  <- apply( tsample, 2, quantile,c(0.025,0.975))


  result  <- list( t= t, variance= vari, lambda = lambda, interval = intval)
# result <- do.call("new", list("bench.sae", result = result) )
} else{
  range  <- diff(interval)
  t      <- mean(interval)
  vari   <- as.vector( (range/4)^2 )
  lambda <- as.vector(1/vari)
  result <- list( t= t, variance= vari, lambda = lambda)
}
return( result )
}






