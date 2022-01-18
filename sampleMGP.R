sampleMGP<-function(theta.jh, delta.h, a1 = 3, a2 = 5){
  # Store the dimensions locally
  p = nrow(theta.jh); K = ncol(theta.jh)
  
  # Sum over the (squared) replicates:
  sum.theta.l = colSums(theta.jh^2) 
  
  # h = 1 case is separate:
  tau.not.1 = cumprod(delta.h)/delta.h[1] 
  delta.h[1] = rgamma(n = 1, shape = a1 + p*K/2, 
                      rate = 1 + 1/2*sum(tau.not.1*sum.theta.l))
  # h > 1:
  if(K > 1){for(h in 2:K){
    tau.not.h = cumprod(delta.h)/delta.h[h]
    delta.h[h] = rgamma(n = 1, shape = a2 + p/2*(K - h + 1), 
                        rate = 1 + 1/2*sum(tau.not.h[h:K]*sum.theta.l[h:K]))
  }}
  delta.h #list(tau.h = cumprod(delta.h), delta.h = delta.h)
}
