#'@title Smooth Binomial data, accounting for nugget effect
#'@param x: a vector of observed number of successes
#'@param size: a vector of number of trials
#'@param sigma: nugget effect
#'@param filter.number,family: wavelet basis, see wavethresh pakcage for more details.
#'@param maxiter: maximum iterations, in general, 1 iteration is enough.
#'@param tol: tolerance to stop iterations.
#'@param approx_p: whether approxmate p by logistic(mu) or using numerical integration to calulcate mean of logit normal distribution.
#'@import logitnorm
#'@import smashr
#'@import ashr
#'@export

smash.gen.binom = function(x,size,sigma=NULL,filter.number = 1,
                           family = "DaubExPhase",approx_p=FALSE){

  if(!ispowerof2(length(x))){
    reflect.x = reflect(x)
    size = reflect(size)$x
    x = reflect.x$x
    idx = reflect.x$idx
  }else{
    idx = 1:length(x)
  }

  n = length(x)
  if(min(x)==0 | max(x==size)==1){
    x_pm = ash(rep(0,n),1,lik=lik_binom(x,size),optmethod='mixSQP',pointmass=T)$result$PosteriorMean
    p_tilde = x/size
    p_tilde[x==0] = x_pm[x==0]
    p_tilde[x==size] = x_pm[x==size]
  }else{
    p_tilde = x/size
  }


  # working data
  st=sqrt(1/(size*p_tilde*(1-p_tilde)))
  y=logit(p_tilde)+(x/size-p_tilde)/(p_tilde*(1-p_tilde))
  sigma.est=sigma




    if(is.null(sigma)){
      fit = NuggetEst(y,st,sigma.est,filter.number = filter.number,family = family)
      sigma.est = sqrt(fit$nugget.est)
    }else{
      fit = smash.gaus(y,sigma=sqrt(st^2+sigma^2),
                       filter.number = filter.number,family = family,post.var = TRUE)
    }

  mu.est = fit$mu.est[idx]
  mu.est.var = fit$mu.est.var[idx]

  if(approx_p){
    p.est = logistic(mu.est)
  }else{
    p.est = apply(rbind(mu.est,mu.est.var),2,function(x) momentsLogitnorm(x[1],sqrt(x[2])))[1,]
  }


  loglik = sum(dnorm(y,mu.est,sd=sqrt(st^2+sigma.est^2),log=TRUE))

  return(list(p.est = p.est,mu.est=mu.est,nugget.est = sigma.est,loglik=loglik))
}

#'@export
logit = function(p){
  log(p/(1-p))
}

#'@export
logistic = function(x){
  1/(1+exp(-x))
}

