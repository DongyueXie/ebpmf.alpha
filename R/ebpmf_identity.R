#'@title Smoothed Poisson Topic Model
#'@description This function fits Poisson Topic Model with smooth Loading or Factors
#'@param X count matrix
#'@param K number of factors/ranks
#'@param S background
#'@param init initialization methods, 'lee','scd' from package NNLM, or 'uniform' randomly initialize; or provide init as a list with L_init and F_init.
#'@param init_loss loss function of the initialization method, either mkl or mse.
#'@param maxiter maximum iterations
#'@param tol stop criteria
#'@param smooth_method splitting or bmsm
#'@param fix_F if TRUE, F will not be updated.
#'@param bmsm_control control parameters of BMSM, see bmsm_control_default()
#'@param ebpm_method point_gamma or two_gamma
#'@param ebpm_control control parameters of ebpm, see ebpm_control_default()
#'@param splitting_control control parameters of smashgen, see splitting_control_default()
#'@param smooth_f,smooth_l whether to get smooth estimate of loadings or factors.
#'@param nugget whether to assume nugget effects
#'@param return_all whether return all outputs or simplified ones
#'@return EL,EF: posterior of loadings and factors
#'@examples
#'set.seed(123)
#'n = 120
#'p = 256
#'K= 3
#'L = matrix(0, nrow=n, ncol=K)
#'FF = matrix(0, nrow=K, ncol=p)
#'L[1:(n/3),1] = 1
#'L[((n/3)+1):(2*n/3),2] = 1
#'L[((2*n/3)+1):n,3] = 1
#'L = L + matrix(runif(n*K,0,0.5),nrow=n)
#'FF[1,1:(p/3)] = 1+10
#'FF[2,((p/3)+1):(2*p/3)] = 1+10
#'FF[3,((2*p/3)+1):p] = 1+10
#'lambda = L %*% FF
#'X = matrix(rpois(n=length(lambda),lambda),nrow=n)
#'image(X)
#'@importFrom NNLM nnmf
#'@import ebpm
#'@import Matrix
#'@import vebpm
#'@importFrom smashrgen pois_smooth_split
#'@export

ebpmf_identity = function(X,K,
                          S=NULL,
                          init = 'scd',
                          init_loss = 'mkl',
                          maxiter=100,
                          tol=1e-8,
                          ebpm.fn=ebpm::ebpm_point_gamma,
                          fix_F = FALSE,
                          printevery=10,
                          verbose=TRUE,
                          convergence_criteria = 'Labs'){


  #browser()
  n = dim(X)[1]
  p = dim(X)[2]
  n_points = n*p

  if(is.null(S)){
    S = 1
  }
  res = init_stm(X,K,init,init_loss)
  EZ = array(dim = c(n,p,K))
  EZ = Calc_EZ(X,K,EZ,res$ql,res$qf)

  #KL = c()
  #KL[1] = mKL(X,tcrossprod(res$ql$El,res$qf$Ef))

  obj = c()
  obj[1] = -Inf

  for(iter in 1:maxiter){
    #b_k_max = 0
    for(k in 1:K){
      res = stm_update_rank1(EZ[,,k],k,ebpm.fn,res,fix_F)
    }
    # Update Z
    EZ = Calc_EZ(X,K,EZ,res$ql,res$qf)

    if(convergence_criteria == 'mKLabs'){
      obj[iter+1] = mKL(X,tcrossprod(res$ql$El,res$qf$Ef))
      if(verbose){
        if(iter%%printevery==0){
          print(sprintf('At iter %d, mKL: %f',iter,obj[iter+1]))
        }
      }
      if(abs(obj[iter+1]-obj[iter])<=tol){
        break
      }

    }

    if(convergence_criteria=='ELBO'){
      obj[iter+1] = calc_stm_obj(X,K,res)
      if(verbose){
        if(iter%%printevery==0){
          print(sprintf('At iter %d, ELBO: %f',iter,obj[iter+1]))
        }
      }
      if((obj[iter+1]-obj[iter])/n_points<tol){
        break
      }
    }

    if(convergence_criteria =='Labs'){
      obj[iter+1] = norm(res$ql$Esmooth_l,type='F')
      if(verbose){
        if(iter%%printevery==0){
          print(sprintf('At iter %d, L norm: %f',iter,obj[iter+1]))
        }
      }
      if(abs(obj[iter+1]-obj[iter])/n_points<tol){
        break
      }
    }


  }
  if(iter==maxiter){
    warning('Reached maximum iterations')
  }

  #lambda_hat = tcrossprod(res$ql$El,res$qf$Ef)
  #lambda_init = L_init%*%F_init

  ldf = poisson_to_multinom(res$qf$Esmooth_f,res$ql$Esmooth_l)
  fit = list(res = res,EL = ldf$L,EF = ldf$FF,d=ldf$s,obj=obj)
  return(fit)

}

calc_stm_obj = function(X,K,res){
  n = nrow(X)
  p = ncol(X)
  val = 0
  qz = calc_qz(X,K,res$ql,res$qf)
  for(k in 1:K){
    val = val + qz[,,k]*(matrix(res$ql$Elogl[,k],nrow=n,ncol=p,byrow=F)+matrix(res$qf$Elogf[,k],nrow=n,ncol=p,byrow=T)-log(qz[,,k]))
  }
  E1 = sum(X*val) - sum(tcrossprod(res$ql$El,res$qf$Ef))

  return(E1+sum(res$Hl)+sum(res$Hf))
  #return(val-sum(tcrossprod(res$ql$El,res$qf$Ef)) +sum(res$Hl)+sum(res$Hf) )
}

calc_H = function(x,s,loglik,pm,pmlog){
  if(is.null(loglik)){
    H = NULL
  }else{
    H = loglik - sum(x*log(s)+x*pmlog-pm*s-lfactorial(x))
  }
  H
}
stm_update_rank1 = function(Z,k,ebpm.fn,res,fix_F){
  if(length(ebpm.fn)==1){
    ebpm.fn.l = ebpm.fn
    ebpm.fn.f = ebpm.fn
  }
  if(length(ebpm.fn)==2){
    ebpm.fn.l = ebpm.fn[[1]]
    ebpm.fn.f = ebpm.fn[[2]]
  }
  # update l
  l_seq = rowSums(Z)
  l_scale = sum(res$qf$Ef[,k])
  fit = ebpm.fn.l(l_seq,l_scale)
  res$ql$El[,k] = fit$posterior$mean
  res$ql$Elogl[,k] = fit$posterior$mean_log
  if(is.null(fit$posterior$mean_smooth)){
    res$ql$Esmooth_l[,k] = fit$posterior$mean
  }else{
    res$ql$Esmooth_l[,k] = fit$posterior$mean_smooth
  }

  res$Hl[k] = calc_H(l_seq,l_scale,fit$log_likelihood,fit$posterior$mean,fit$posterior$mean_log)
  #res$Esmooth_l = lk_hat$Esmooth
  if(!fix_F){
    # update f
    f_seq = colSums(Z)
    f_scale = sum(res$ql$El[,k])
    fit = ebpm.fn.f(f_seq,f_scale)
    res$qf$Ef[,k] = fit$posterior$mean
    res$qf$Elogf[,k] = fit$posterior$mean_log
    if(is.null(fit$posterior$mean_smooth)){
      res$qf$Esmooth_f[,k] = fit$posterior$mean
    }else{
      res$qf$Esmooth_f[,k] = fit$posterior$mean_smooth
    }
    res$Hf[k] = calc_H(f_seq,f_scale,fit$log_likelihood,fit$posterior$mean,fit$posterior$mean_log)

  }

  return(res)

}

#'@title initialize the stm model
#'@param X input data matrix
#'@param K number of topics
#'@param init init methods, or a list of init L and F
#'@param init_loss mkl or mse
#'@importFrom NNLM nnmf
#'@export
init_stm = function(X,K,init,init_loss,maxiter_init = 50){

  n = nrow(X)
  if(is.list(init)){
    L_init = init$L_init
    F_init = init$F_init

    if(is.null(L_init)){
      X_init_fit = nnmf(as.matrix(X),K,method='lee',
                              loss='mse',show.warning = F,
                              init = list(H=t(F_init)),
                              verbose = F,max.iter = maxiter_init)
      L_init = X_init_fit$W
    }

  }else{

    if(init%in%c('scd','lee')){
      X_init_fit = NNLM::nnmf(as.matrix(X),K,method=init,loss=init_loss,show.warning = F,verbose = F,max.iter = maxiter_init)
      L_init = X_init_fit$W
      F_init = t(X_init_fit$H)
    }
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
  }

  # adjust scale of L and F, mainly for stability.
  ratio = adjLF(L_init,F_init)
  L_init = ratio$L_init
  F_init = ratio$F_init

  gl = list()
  gf = list()

  ql = list(El = L_init, Elogl = log(L_init+1e-10), Esmooth_l=L_init)
  qf = list(Ef = F_init, Elogf = log(F_init+1e-10), Esmooth_f=F_init)

  return(list(ql=ql,
              qf=qf,
              gl=gl,
              gf=gf,
              Hl = rep(0,K),
              Hf = rep(0,K)))

}


#' #'@title Default parameters of ebpm
#' #'@export
#' ebpm_control_default = function(){
#'   list(pi0 = 'estimate',
#'        g_init = NULL,
#'        fix_g = FALSE,
#'        control =  NULL)
#' }

#'
#' #'@title Default parameters of smash gen
#' #'@param filter.number,family wavelet basis, see wavethresh pakcage for more details.
#' #'@export
#' splitting_control_default = function(){
#'   list(Eb_init = NULL,
#'        sigma2_init = NULL,
#'        est_sigma2 = TRUE,
#'        maxiter = 100,
#'        tol=1e-5,
#'        filter.number = 1,
#'        family = 'DaubExPhase',
#'        verbose=FALSE,
#'        printevery = 10,
#'        ebnm_params=list(mode=0),
#'        optim_method='L-BFGS-B')
#' }
#'



Calc_EZ = function(X,K,EZ,ql_hat,qf_hat){
  n = nrow(X)
  p = ncol(X)
  for(k in 1:K){
    EZ[,,k] = outer(ql_hat$Elogl[,k], qf_hat$Elogf[,k], "+")
  }
  EZ = softmax3d(EZ)
  EZ = as.vector(EZ)*as.vector(X)
  dim(EZ) = c(n,p,K)
  EZ
}

calc_qz = function(X,K,ql,qf){
  n = nrow(X)
  p = ncol(X)
  qz = array(dim = c(n,p,K))
  for(k in 1:K){
    qz[,,k] = outer(ql$Elogl[,k], qf$Elogf[,k], "+")
  }
  return(softmax3d(qz))
}


