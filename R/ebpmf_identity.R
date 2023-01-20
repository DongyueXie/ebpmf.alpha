#'@title Smoothed Poisson Topic Model
#'@description This function fits Poisson Topic Model with smooth Loading or Factors
#'@param X count matrix
#'@param K number of factors/ranks
#'@param init initialization methods, 'lee','scd' from package NNLM, or 'uniform' randomly initialize; or provide init as a list with L_init and F_init.
#'@param init_loss loss function of the initialization method, either mkl or mse.
#'@param maxiter maximum iterations
#'@param tol stop criteria
#'@param ebpm.fn specify functions to use for solving the poisson subproblems
#'@param fix_F if TRUE, F will not be updated.
#'@param smooth_l,smooth_f whether smooth l or f, must match the functions in ebpm.fn
#'@param warm_start whether warm starts the ebpm. Useful for when smooth_l or f is true.
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
                          init = 'fasttopics',
                          maxiter=100,
                          maxiter_init = 30,
                          tol=1e-3,
                          ebpm.fn=c(ebpm::ebpm_point_gamma,smashrgen::pois_smooth_split),
                          fix_F = FALSE,
                          smooth_l = FALSE,
                          smooth_f = TRUE,
                          smooth_control=list(),
                          printevery=10,
                          verbose=TRUE,
                          convergence_criteria = 'mKLabs'){


  start_time = Sys.time()
  #browser()
  n = dim(X)[1]
  p = dim(X)[2]
  n_points = n*p

  if(verbose){
    cat('initializing...')
    cat('\n')
  }
  res = ebpmf_identity_init(X,K,init,maxiter_init)
  EZ = array(dim = c(n,p,K))
  EZ = Calc_EZ(X,K,EZ,res$ql,res$qf)

  #KL = c()
  #KL[1] = mKL(X,tcrossprod(res$ql$El,res$qf$Ef))

  obj = c()
  obj[1] = -Inf

  smooth_control = modifyList(smooth_control_default(),smooth_control,keep.null = TRUE)

  if(verbose){
    cat('running iterations')
    cat('\n')
  }
  for(iter in 1:maxiter){
    #b_k_max = 0
    for(k in 1:K){
      res = stm_update_rank1(EZ[,,k],k,ebpm.fn,res,fix_F,smooth_l,smooth_f,smooth_control)
    }
    # Update Z
    EZ = Calc_EZ(X,K,EZ,res$ql,res$qf)

    if(convergence_criteria == 'mKLabs'){
      obj[iter+1] = mKL(X,tcrossprod(res$ql$Esmooth_l,res$qf$Esmooth_f))
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
    message('Reached maximum iterations')
  }

  #lambda_hat = tcrossprod(res$ql$El,res$qf$Ef)
  #lambda_init = L_init%*%F_init

  ldf = poisson_to_multinom(res$qf$Esmooth_f,res$ql$Esmooth_l)
  fit = list(res = res,EL = ldf$L,EF = ldf$FF,d=ldf$s,obj=obj,run_time = difftime(Sys.time(),start_time,units='auto'))
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
    H = 0
  }else{
    H = loglik - sum(x*log(s)+x*pmlog-pm*s-lfactorial(x))
  }
  H
}
stm_update_rank1 = function(Z,k,ebpm.fn,res,fix_F,smooth_l,smooth_f,smooth_control){
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
    if(smooth_f){
      #print(res$gf)
      fit = ebpm.fn.f(f_seq,f_scale,
                      m_init=smooth_control$m_init,
                      sigma2_init=res$gf$sigma2[k],
                      #smooth_init = res$qf$Esmooth_f[,k],
                      maxiter=smooth_control$maxiter,
                      wave_trans=smooth_control$wave_trans,
                      ndwt_method=smooth_control$ndwt_method,
                      make_power_of_2=smooth_control$make_power_of_2,
                      ash_pm_init_for0 = smooth_control$ash_pm_init_for0,
                      vga_tol=smooth_control$vga_tol,
                      tol = smooth_control$tol,
                      warmstart=smooth_control$warm_start)
    }else{
      fit = ebpm.fn.f(f_seq,f_scale)
    }
    res$qf$Ef[,k] = fit$posterior$mean
    res$qf$Elogf[,k] = fit$posterior$mean_log
    #print(fit$fitted_g$sigma2)
    res$gf$sigma2[k] = fit$fitted_g$sigma2
    if(is.null(fit$posterior$mean_smooth)){
      res$qf$Esmooth_f[,k] = fit$posterior$mean
    }else{
      res$qf$Esmooth_f[,k] = fit$posterior$mean_smooth
    }
    res$Hf[k] = calc_H(f_seq,f_scale,fit$log_likelihood,fit$posterior$mean,fit$posterior$mean_log)

  }

  return(res)

}




#' #'@title Default parameters of ebpm
#' #'@export
#' ebpm_control_default = function(){
#'   list(pi0 = 'estimate',
#'        g_init = NULL,
#'        fix_g = FALSE,
#'        control =  NULL)
#' }


#'@title Default parameters of smooth split
smooth_control_default = function(){
  list(wave_trans='ndwt',
       ndwt_method = "ti.thresh",
       m_init='vga',
       maxiter=30,
       make_power_of_2='extend',
       ash_pm_init_for0 = FALSE,
       vga_tol=1e-3,
       tol = 1e-2,
       warm_start=TRUE)
}




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


