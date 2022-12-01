
#'@title simulate Poisson matrix
#'@details The DGP is
#'\deqn{Y\sim Poisson(\exp(l_0f_0'+LF'+E)),}
#'\deqn{l_{0i}\sim Uniform(),f_{0j}\sim Uniform(),}
#'\deqn{l_{ki}\sim N(0,d_k),f_{kj}\sim \pi_0\delta_0 + (1-\pi_0)N(0,\sigma^2_1),}
#'\deqn{e_{ij}\sim N(0,\sigma^2).}
#'@export
sim_data_log = function(n,p,K=3,d=NULL,
                        add_background = TRUE,
                        background_unif_range = c(1,2),
                        pi0=0.8,var1 = 1,
                        var_e = 0,
                        n_simu = 10,
                        seed=12345){
  set.seed(seed)
  if(is.null(d)){
    d = rep(1,K)
  }

  Y = array(dim=c(n,p,n_simu))
  Loading = array(dim=c(n,K,n_simu))
  Factor = array(dim=c(p,K,n_simu))
  L0 = matrix(nrow=n,ncol=n_simu)
  F0 = matrix(nrow=p,ncol=n_simu)

  for(i in 1:n_simu){
    # f is sparse
    FF = matrix(0,nrow=p,ncol=K)
    for(k in 1:K){
      FF[,k] = rnorm(p,0,sd=sqrt(var1))*rbinom(p,1,1-pi0)
    }
    # L controls the norms
    L = matrix(0,nrow=n,ncol=K)
    for(k in 1:K){
      L[,k] = rnorm(n,0,sqrt(d[k]))
    }
    if(add_background){
      l0 = runif(n,min=background_unif_range[1],background_unif_range[2])
      f0 = runif(p,min=background_unif_range[1],background_unif_range[2])
      L0[,i] = l0
      F0[,i] = f0
      S0 = tcrossprod(l0,f0)
    }else{
      S0 = 1
    }
    Factor[,,i] = FF
    Loading[,,i] = L
    Lambda = S0*exp(tcrossprod(L,FF)+matrix(rnorm(n*p,0,sqrt(var_e)),nrow=n,ncol=p))
    Y[,,i] = matrix(rpois(n*p,Lambda),nrow=n,ncol=p)

  }

  return(list(Y=Y,Factor=Factor,Loading=Loading,L0=L0,F0=F0,params = list(n=n,p=p,n_simu=n_simu,
                                                              d=d,
                                                              add_background = add_background,
                                                              background_unif_range = background_unif_range,
                                                              pi0=pi0,var1 = var1,
                                                              var_e = var_e)))

}

#'@title simulate simple Poisson matrix
#'@export
sim_data_log_simple = function(N,p){
  K = 2
  Ftrue = matrix(0,nrow=p,ncol=K)
  Ftrue[1:20,1] = 1
  Ftrue[21:40,2] = 1
  Ltrue = matrix(rnorm(N*K), ncol=K)

  l0 = runif(N,1,2)
  f0 = runif(p,1,2)
  S0 = tcrossprod(l0,f0)

  Lambda = S0*exp(tcrossprod(Ltrue,Ftrue))

  Y = matrix(rpois(N*p,Lambda),nrow=N,ncol=p)
  return(list(Y=Y,l0=l0,f0=f0,Ftrue = Ftrue,Ltrue=Ltrue))
}




