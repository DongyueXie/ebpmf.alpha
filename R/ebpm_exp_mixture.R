#'@title Empirical Bayes Poisson Mean
#'@description Modified from Zihao Wang's ebpm package
#'@export
#'


ebpm_exp_mixture <- function(x,s = 1, scale = c(1e-10,1e-5,0.01,0.1,1,10,100,1e5), point_mass=F,
                             nullweight=1000, weight = rep(1,length(x)),
                             g_init = NULL, fix_g = FALSE,
                             m = 2, control =  NULL, low = NULL,d=NULL,shape=1){
  n=length(x)
  if(length(s) == 1){s = replicate(length(x),s)}
  if(is.null(control)){control = mixsqp_control_defaults()}
  if(is.null(g_init)){
    fix_g = FALSE ## then automatically unfix g if specified so
    if(identical(scale, "estimate")){scale <- select_grid_exponential(x,s,m,d, low,shape)
    }else{
      scale = list(scale=scale,shape=rep(shape,length(scale)))
    }
    g_init = scale2gammamix_init(scale,point_mass)
  }

  if(!fix_g){ ## need to estimate g_init
    b = 1/g_init$scale ##  from here use gamma(shape = a, rate = b)  where E = a/b
    a = g_init$shape
    tmp <-  compute_L(x,s,a, b,point_mass)
    L =  tmp$L
    l_rowmax = tmp$l_rowmax
    if(point_mass){x0 = c(g_init$pi0,g_init$pi)}else{x0 = g_init$pi}
    if(!is.null(nullweight)){
      Lnull = rbind(c(1,rep(0,ncol(L)-1)),L)
      weight = c(nullweight-1,weight)
      fit <- mixsqp(Lnull, weight,x0 = x0, control = control)
    }else{
      fit <- mixsqp(L, weight,x0 = x0, control = control)
    }
    pi = fit$x
    pi = pi/sum(pi) ## seems that some times pi does not sum to one
  }
  else{
    if(point_mass){
      pi = c(g_init$pi0,g_init$pi)
    }else{
      pi = g_init$pi
    }
    a = g_init$shape
    b = 1/g_init$scale
    ## compute loglikelihood
    tmp <-  compute_L(x,s,a, b,point_mass)
    L =  tmp$L
    l_rowmax = tmp$l_rowmax
  }
  fitted_g = gammamix(pi = pi, shape = a,  scale  = 1/b,point_mass)

  log_likelihood = sum(log(exp(l_rowmax) * L %*%  pi))

  cpm = outer(x,a,  "+")/outer(s, b, "+")
  if(point_mass){cpm = cbind(rep(0,n),cpm)}
  Pi_tilde = t(t(L) * pi)
  Pi_tilde = Pi_tilde/rowSums(Pi_tilde)
  lam_pm = rowSums(Pi_tilde * cpm)

  c_log_pm = digamma(outer(x,a,  "+")) - log(outer(s, b, "+"))
  if(point_mass){
    lam_log_pm = rowSums(Pi_tilde[,-1] * c_log_pm)
    lam_log_pm[x==0] = -Inf
  }else{
    lam_log_pm = rowSums(Pi_tilde * c_log_pm)
  }
  posterior = data.frame(mean = lam_pm, mean_log = lam_log_pm)
  return(list(fitted_g = fitted_g,
              posterior = posterior,
              log_likelihood = log_likelihood,
              Pi_tilde=Pi_tilde,
              tmp=tmp))
}

geom_seq <- function(low, up, m){
  N =  ceiling(log(up/low)/log(m)) + 1
  out  = low*m^(seq(1,N, by = 1)-1)
  return(out)
}


prune_fitted_g_ebpm = function(fitted_g,thresh=1e-10){
  rm_idx = which(fitted_g$pi<thresh)
  fitted_g$pi = fitted_g$pi[-rm_idx]
  fitted_g$shape = fitted_g$shape[-rm_idx]
  fitted_g$scale = fitted_g$scale[-rm_idx]
  fitted_g
}

## select grid for b_k
select_grid_exponential <- function(x, s, m = 2, d = NULL, low = NULL,shape){
  ## mu_grid: mu =  1/b is the exponential mean
  xprime = x
  xprime[x == 0] = xprime[x == 0] + 1
  mu_grid_min =  0.05*min(xprime/s)
  mu_grid_max = 2*max(x/s)
  if(is.null(m)){
    if(is.null(d)){m = 2}
    else{m = ceiling((mu_grid_max/mu_grid_min)^(1/(d-1)))}
  }
  if(!is.null(low)){mu_grid_min = min(low, mu_grid_min)}
  if(mu_grid_min<1e-10){mu_grid_min=1e-10}
  mu_grid = geom_seq(mu_grid_min, mu_grid_max, m)
  a = rep(shape, length(mu_grid))
  return(list(shape = a, scale = mu_grid))
}



#' @export get_uniform_mixture
get_uniform_mixture <- function(x, s, grid_res = NULL, m = 2, low = NULL,point_mass){
  if(is.null(grid_res)){
    grid_res = select_grid_exponential(x = x, s = s, m = m, low = low)
  }
  shape = grid_res$shape
  scale = grid_res$scale
  n = length(shape)
  pi = replicate(n, 1/n)
  g = gammamix(pi = pi, shape = shape, scale = scale,point_mass)
  return(g)
}

## compute L matrix from data and selected grid
## L_ik = NB(x_i; a_k, b_k/b_k + s_i)
## but for computation in mixsqr, we can simplyfy it for numerical stability
compute_L <- function(x, s, a, b,point_mass){
  prob = 1 - s/outer(s,b, "+")
  l = dnbinom_cts_log(x,a,prob = prob) ##
  l_rowmax  = apply(l,1,max)
  if(point_mass){
    l0 = cbind(log(c(x==0)),l)
    L = exp(l0 -  l_rowmax)
  }else{
    L = exp(l -  l_rowmax)
  }

  return(list(L = L, l_rowmax = l_rowmax))
}


# it is equivalent to dnbinom in R wiht log = T when X is integer; I allow  it  to compute when x is not integer
dnbinom_cts_log <- function(x, a, prob){
  tmp = x*log(1-prob)
  tmp[x == 0] = 0 ## R says 0*-Inf = NaN
  out = t(t(log(prob)) * a) + tmp + lgamma(outer(x, a, "+")) - lgamma(x+1)
  out = t(t(out) - lgamma(a))
  return(out)
}

gammamix <- function(pi, shape, scale,point_mass) {
  if(point_mass){
    structure(list(pi = pi[-1], shape = shape, scale = scale, pi0 = pi[1]), class="gammamix")
  }else{
    structure(list(pi = pi, shape = shape, scale = scale), class="gammamix")
  }
}

scale2gammamix_init <- function(scale,point_mass){
  n = length(scale$shape) + point_mass
  pi_init = replicate(n, 1)/n
  return(gammamix(pi = pi_init, shape = scale$shape, scale =  scale$scale,point_mass))
}

mixsqp_control_defaults <- function() {
  return(list(verbose = F))
}
