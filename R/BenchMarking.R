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
















dim( THETA )
theta_b = expit( apply(THETA,1,mean) )
sum(theta_b)
t = 3
w = rbinom(length(theta_b),10,0.5) +1
w = w/sum(w)
phi=rbinom(length(theta_b),10,0.5) +1
phi = phi/sum(phi)





sum( theta_b*w )
BM(theta_b,2.5,w,phi)
BM(theta_b,2.5,w,phi,1000)
BM(theta_b,2.5,w,phi,100)

sum( BM(theta_b,2.5,w,phi)*w )
sum( BM(theta_b,2.5,w,phi,1000)*w )
sum( BM(theta_b,2.5,w,phi,100)*w )
2.5





