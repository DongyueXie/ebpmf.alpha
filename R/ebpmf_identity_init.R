#'@title initialize the ebpmf identity model
#'@param X input data matrix
#'@param K number of topics
#'@param init init methods, or a list of init L and F
#'@importFrom fastTopics fit_poisson_nmf
#'@export
ebpmf_identity_init = function(X,K,init,maxiter_init = 50){

  n = nrow(X)
  if(is.list(init)){
    L_init = init$L_init
    F_init = init$F_init

    if(is.null(L_init)&!is.null(F_init)){
      init_fasttopic = init_poisson_nmf(X,F=F_init,init.method='random')
      init_fit = fit_poisson_nmf(X,numiter = maxiter_init,fit0 = init_fasttopic)
      L_init = init_fit$L
    }

  }else{

    if(init == 'uniform'){
      L_init = matrix(runif(n*K),nrow=n,ncol=K)
      F_init = matrix(runif(K*p),nrow=p,ncol=K)
      ratio = median(X)/(median(L_init)*median(F_init))
      L_init = L_init*sqrt(ratio)
      F_init = F_init*sqrt(ratio)
    }
    if(init == 'kmeans'){
      kmeans.init=kmeans(as.matrix(X),K,nstart=5)
      L_init = rep(1,n)%o%normalize(as.vector(table(kmeans.init$cluster)))
      F_init = t(kmeans.init$centers)
      row.names(F_init)=NULL
    }
    if(init == 'fasttopics'){
      init_fasttopic = fit_poisson_nmf(X,K,numiter = maxiter_init)
      L_init = init_fasttopic$L
      F_init = init_fasttopic$F
    }
  }

  # adjust scale of L and F, mainly for stability.
  ratio = adjLF(L_init,F_init)
  L_init = ratio$L_init
  F_init = ratio$F_init

  gl = list(sigma2=NULL)
  gf = list(sigma2=NULL)

  ql = list(El = L_init, Elogl = log(L_init+1e-10), Esmooth_l=L_init)
  qf = list(Ef = F_init, Elogf = log(F_init+1e-10), Esmooth_f=F_init)

  return(list(ql=ql,
              qf=qf,
              gl=gl,
              gf=gf,
              Hl = rep(0,K),
              Hf = rep(0,K)))

}




#'
#' #'@title initialize the ebpmf identity model
#' #'@param X input data matrix
#' #'@param K number of topics
#' #'@param init init methods, or a list of init L and F
#' #'@param init_loss mkl or mse
#' #'@importFrom NNLM nnmf
#' #'@export
#' ebpmf_identity_init = function(X,K,init,init_loss,maxiter_init = 50){
#'
#'   n = nrow(X)
#'   if(is.list(init)){
#'     L_init = init$L_init
#'     F_init = init$F_init
#'
#'     if(is.null(L_init)){
#'       X_init_fit = nnmf(as.matrix(X),K,method='lee',
#'                         loss='mse',show.warning = F,
#'                         init = list(H=t(F_init)),
#'                         verbose = F,max.iter = maxiter_init)
#'       L_init = X_init_fit$W
#'     }
#'
#'   }else{
#'
#'     if(init%in%c('scd','lee')){
#'       X_init_fit = NNLM::nnmf(as.matrix(X),K,method=init,loss=init_loss,show.warning = F,verbose = F,max.iter = maxiter_init)
#'       L_init = X_init_fit$W
#'       F_init = t(X_init_fit$H)
#'     }
#'     if(init == 'uniform'){
#'       L_init = matrix(runif(n*K),nrow=n,ncol=K)
#'       F_init = matrix(runif(K*p),nrow=p,ncol=K)
#'       ratio = median(X)/(median(L_init)*median(F_init))
#'       L_init = L_init*sqrt(ratio)
#'       F_init = F_init*sqrt(ratio)
#'     }
#'     if(init == 'kmeans'){
#'       kmeans.init=kmeans(as.matrix(X),K,nstart=5)
#'       L_init = rep(1,n)%o%normalize(as.vector(table(kmeans.init$cluster)))
#'       F_init = t(kmeans.init$centers)
#'       row.names(F_init)=NULL
#'     }
#'   }
#'
#'   # adjust scale of L and F, mainly for stability.
#'   ratio = adjLF(L_init,F_init)
#'   L_init = ratio$L_init
#'   F_init = ratio$F_init
#'
#'   gl = list(sigma2=NULL)
#'   gf = list(sigma2=NULL)
#'
#'   ql = list(El = L_init, Elogl = log(L_init+1e-10), Esmooth_l=L_init)
#'   qf = list(Ef = F_init, Elogf = log(F_init+1e-10), Esmooth_f=F_init)
#'
#'   return(list(ql=ql,
#'               qf=qf,
#'               gl=gl,
#'               gf=gf,
#'               Hl = rep(0,K),
#'               Hf = rep(0,K)))
#'
#' }
#'
