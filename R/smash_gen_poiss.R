#'@title Smooth Poisson data, accounting for nugget effect
#'@param x: observed Poisson sequence
#'@param sigma: nugget effect
#'@param s: Scale factor for Poisson observations: y~Pois(scale*lambda).
#'@param method: smoothiing method, smash or ti.thresh. When n is large, ti.thresh is much faster.
#'@param ash.pm: whehter use ash posterior mean approxiamtion if x=0. If not x = x+eps
#'@param eps: if x=0, x = x + eps
#'@param filter.number,family: wavelet basis, see wavethresh pakcage for more details.
#'@param maxiter: maximum iterations, in general, 1 iteration is enough.
#'@param tol: tolerance to stop iterations.
#'@import smashr
#'@import ashr
#'@export

smash.gen.poiss = function(x,sigma=NULL,s=1,method='ti.thresh',ash.pm=FALSE,eps=0.01, filter.number = 1,family = "DaubExPhase"){

  if(!ispowerof2(length(x))){
    reflect.x = reflect(x)
    x = reflect.x$x
    idx = reflect.x$idx
  }else{
    idx = 1:length(x)
  }

  n = length(x)
  if(min(x)<1){
    if(ash.pm){
      x_pm = ash(rep(0,n),1,lik=lik_pois(x,scale=s),optmethod='mixSQP',pointmass=T)$result$PosteriorMean
      lambda_tilde = x/s
      lambda_tilde[x<1] = x_pm[x<1]
    }else{
      lambda_tilde = (x+eps)/s
    }
  }else{
    lambda_tilde = x/s
  }


  # working data
  st=sqrt(1/(s*lambda_tilde))
  y=log(lambda_tilde)+(x-s*lambda_tilde)/(s*lambda_tilde)



  #estimate nugget effect and estimate mu

  sigma.est = sigma



    if(is.null(sigma)){
      fit = NuggetEst(y,st,sigma.est,method,filter.number = filter.number,family = family)
      sigma.est = sqrt(fit$nugget.est)
      mu.est = fit$mu.est[idx]
    }else{
      if(method=='smash'){
        fit = smash.gaus(y,sigma=sqrt(st^2+sigma^2),
                         filter.number = filter.number,family = family,
                         post.var = TRUE,return.loglr = TRUE)
        mu.est = fit$mu.est[idx]
      }
      if(method == 'ti.thresh'){
        fit = ti.thresh(y,sigma=sqrt(st^2+sigma^2),filter.number = filter.number,family = family)
        mu.est = fit[idx]
      }
    }



  if(method=='smash'){
    mu.est.var = fit$mu.est.var[idx]
    lambda.est = exp(mu.est+mu.est.var/2)
  }
  if(method=='ti.thresh'){
    lambda.est = exp(mu.est)
  }


  #lambda.est = exp(mu.est)

  #loglik = fit$loglik


  return(list(lambda.est=lambda.est,mu.est=mu.est,nugget.est=sigma.est))
}



normaleqn=function(nug,y,mu,st){
  return(sum((y-mu)^2/(nug+st^2)^2)-sum(1/(nug+st^2)))
}

NuggetEst=function(y,st,nug.init=NULL,method,filter.number,family){
  #initialize nugget effect sigma^2
  n=length(y)
  if(is.null(nug.init)){
    x.m=c(y[n],y,y[1])
    st.m=c(st[n],st,st[1])
    nug.init = ((x.m[2:n]-x.m[3:(n+1)])^2+(x.m[2:n]-x.m[1:(n-1)])^2-2*st.m[2:n]^2-st.m[1:(n-1)]^2-st.m[3:(n+1)]^2)/4
    nug.init = nug.init[nug.init>0&nug.init<var(y)]
    nug.init = mean(nug.init)
  }
  #given st and nug to estimate mean\

  if(method == 'smash'){
    mean.est=smash.gaus(y,sigma=sqrt(st^2+nug.init),filter.number = filter.number,family = family)
    #given mean estimate nugget effect
    nug.est=uniroot(normaleqn,c(-1e6,1e6),y=y,mu=mean.est,st=st)$root
    nug.est = max(c(0,nug.est))
    mean.est=smash.gaus(y,sigma=sqrt(st^2+nug.est),
                        filter.number = filter.number,family = family,
                        post.var = TRUE,return.loglr = TRUE)
    #if wanna mean estimation output, then estiamte mean again
    return(list(mu.est=mean.est$mu.est,mu.est.var=mean.est$mu.est.var,
                nugget.est=nug.est,loglik=mean.est$loglik))
  }else{
    mean.est = ti.thresh(y,sigma=sqrt(st^2+nug.init),filter.number=filter.number,family=family)
    nug.est=uniroot(normaleqn,c(-1e6,1e6),y=y,mu=mean.est,st=st)$root
    nug.est = max(c(0,nug.est))
    mean.est = ti.thresh(y,sigma=sqrt(st^2+nug.init),filter.number=filter.number,family=family)
    return(list(mu.est = mean.est,nugget.est=nug.est))
  }



}


ispowerof2 <- function (x){
  x >= 1 & 2^ceiling(log2(x)) == x
}




