## Benchmarking function

BM = function(theta_b,t,w,phi,lambda=NA){
r = w/phi
s = sum(w*r)
if( is.na(lambda) ){
theta_BM = theta_b + ( (s )^(-1) )*(t - sum(w*theta_b) )*r
} else{
theta_BM = theta_b + ( (s + lambda^(-1))^(-1) )*(t - sum(w*theta_b) )*r}
return(theta_BM)
}
