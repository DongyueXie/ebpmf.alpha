#'@title Smoothed Poisson Topic Model
#'@description This function fits Poisson Topic Model with smooth Loading or Factors
#'@param X count matrix
#'@param K number of factors/ranks
#'@param init initialization methods, 'fasttopics' or 'uniform' randomly initialize; or provide init as a list with L_init and F_init.
#'@param maxiter,maxiter_init maximum iterations
#'@param tol stop criteria
#'@param ebpm.fn specify functions to use for solving the poisson subproblems
#'@param fix_F if TRUE, F will not be updated.
#'@param smooth_F whether smooth l or f, must match the functions in ebpm.fn
#'@param smooth_control a list. ebpmf_identity_smooth_control_default() gives default settings.
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
#'@import ebpm
#'@import Matrix
#'@import vebpm
#'@importFrom smashrgen ebps
#'@importFrom Rfast rowsums
#'@export

ebpmf_identity = function(X,K,
                          init = 'fasttopics',
                          maxiter=50,
                          maxiter_init = 30,
                          tol=1e-3,
                          ebpm.fn=c(ebpm::ebpm_point_gamma,smashrgen::ebps),
                          fix_F = FALSE,
                          smooth_F = TRUE,
                          smooth_control=list(),
                          printevery=10,
                          verbose=TRUE,
                          convergence_criteria = 'mKLabs'){


  start_time = Sys.time()
  #browser()
  n = dim(X)[1]
  p = dim(X)[2]
  n_points = n*p

  X = Matrix(X,sparse = T)
  x = summary(X)
  non0_idx = cbind(x$i,x$j)

  smooth_control = modifyList(ebpmf_identity_smooth_control_default(),smooth_control,keep.null = TRUE)
  if(length(ebpm.fn)==1){
    ebpm.fn.l = ebpm.fn
    ebpm.fn.f = ebpm.fn
  }
  if(length(ebpm.fn)==2){
    ebpm.fn.l = ebpm.fn[[1]]
    ebpm.fn.f = ebpm.fn[[2]]
  }

  if(verbose){
    cat('initializing loadings and factors...')
    cat('\n')
  }

  res = ebpmf_identity_init(X,K,init,maxiter_init)
  alpha = res$ql$Elogl[x$i,] + res$qf$Elogf[x$j,]
  exp_offset = rowMaxs(alpha)
  alpha = alpha - outer(exp_offset,rep(1,K),FUN='*')
  alpha = exp(alpha)
  alpha = alpha/rowsums(alpha)

  # init for smooth curve

  if(smooth_F){
    if(verbose){
      cat('initializing smooth factors...')
      cat('\n')
    }
    res = ebpmf_identity_init_smooth(res,K,p,x,alpha,smooth_control,maxiter_init,ebpm.fn.f)
  }


  obj = c()
  obj[1] = -Inf

  if(verbose){
    cat('running iterations')
    cat('\n')
  }
  if(smooth_F){
    res$ldf = poisson_to_multinom(res$qf$Ef_smooth,res$ql$El)
  }else{
    res$ldf = poisson_to_multinom(res$qf$Ef,res$ql$El)
  }
  for(iter in 1:maxiter){

    for(k in 1:K){
      Ez = calc_EZ(x, alpha[,k])
      res = stm_update_rank1(Ez$rs,Ez$cs,k,ebpm.fn.l,ebpm.fn.f,res,fix_F,smooth_F,smooth_control)
    }
    # Update Z
    # EZ = Calc_EZ(X,K,EZ,res$ql,res$qf)
    alpha = res$ql$Elogl[x$i,] + res$qf$Elogf[x$j,]
    exp_offset = rowMaxs(alpha)
    alpha = alpha - outer(exp_offset,rep(1,K),FUN='*')
    alpha = exp(alpha)
    alpha = alpha/rowsums(alpha)

    if(convergence_criteria == 'mKLabs'){
      if(smooth_F){
        obj[iter+1] = mKL(x$x,tcrossprod(res$ql$El,res$qf$Ef_smooth)[non0_idx])
      }else{
        obj[iter+1] = mKL(x$x,tcrossprod(res$ql$El,res$qf$Ef)[non0_idx])
      }

      if(verbose){
        if(iter%%printevery==0){
          cat(sprintf('At iter %d, mKL(X,LF) = %f',iter,obj[iter+1]))
          cat('\n')
        }
      }
      if(abs(obj[iter+1]-obj[iter])<=tol){
        break
      }
    }

    if(convergence_criteria=='ELBO'){
      obj[iter+1] = calc_stm_obj(x,n,p,K,res,non0_idx)
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
      if(smooth_F){
        ldf = poisson_to_multinom(res$qf$Ef_smooth,res$ql$El)
      }else{
        ldf = poisson_to_multinom(res$qf$Ef,res$ql$El)
      }
      obj[iter+1] = norm(res$ldf$L - ldf$L,type='F')
      res$ldf = ldf
      if(verbose){
        if(iter%%printevery==0){
          print(sprintf('At iter %d, ||L new - L|| F norm: %f',iter,obj[iter+1]))
        }
      }
      if(obj[iter+1]/n_points<tol){
        break
      }
    }
  }
  if(iter==maxiter){
    message('Reached maximum iterations')
  }

  # calc elbo(approximated)
  if(verbose){
    cat('wrapping-up')
    cat('\n')
  }
  if(smooth_F & !fix_F){
    res = calc_approx_elbo_F(x,alpha,K,ebpm.fn.f,res,smooth_control)
  }
  elbo = calc_stm_obj(x,n,p,K,res,non0_idx)


  if(smooth_F){
    ldf = poisson_to_multinom(res$qf$Ef_smooth,res$ql$El)
  }else{
    ldf = poisson_to_multinom(res$qf$Ef,res$ql$El)
  }
  fit = list(EL = ldf$L,
             EF = ldf$FF,
             elbo=elbo,
             d=ldf$s,
             obj=obj,
             res = res,
             run_time = difftime(Sys.time(),start_time,units='auto'))
  return(fit)
}

# when an approxiamted elbo is needed. This is useful when smooth F.
calc_approx_elbo_F = function(x,alpha,K,ebpm.fn.f,res,ebps_control){

  for(k in 1:K){
    Ez = calc_EZ(x, alpha[,k])
    fit = ebpm.fn.f(Ez$cs,sum(res$ql$El[,k]),
                    g_init = list(sigma2 = res$gf$sigma2[k]),
                    q_init = list(m=res$qf$Elogf[,k],smooth = res$qf$Elogf_smooth[,k]),
                    general_control = list(maxiter=1,
                                           maxiter_vga = ebps_control$maxiter_vga,
                                           make_power_of_2=ebps_control$make_power_of_2,
                                           vga_tol=ebps_control$vga_tol,
                                           tol = ebps_control$tol),
                    smooth_control = list(wave_trans='dwt',
                                          ndwt_method = ebps_control$ndwt_method,
                                          filter.number = ebps_control$filter.number,
                                          family = ebps_control$family,
                                          ebnm_params=ebps_control$ebnm_params,
                                          warmstart=ebps_control$warmstart))
    res$Hf[k] = calc_H(Ez$cs,sum(res$ql$El[,k]),fit$log_likelihood,fit$posterior$mean,fit$posterior$mean_log)
  }

  res

}

calc_stm_obj = function(x,n,p,K,res,non0_idx){
  val = 0
  qz = calc_qz(n,p,K,res$ql,res$qf)
  for(k in 1:K){
    val = val + qz[,,k]*(matrix(res$ql$Elogl[,k],nrow=n,ncol=p,byrow=F)+matrix(res$qf$Elogf[,k],nrow=n,ncol=p,byrow=T)-log(qz[,,k]))
  }
  E1 = sum(x$x*val[non0_idx]) - sum(tcrossprod(res$ql$El,res$qf$Ef))

  return(E1+sum(res$Hl)+sum(res$Hf))
}

calc_H = function(x,s,loglik,pm,pmlog){
  if(is.null(loglik)){
    H = 0
  }else{
    H = loglik - sum(x*log(s)+x*pmlog-pm*s-lfactorial(x))
  }
  H
}

stm_update_rank1 = function(l_seq,f_seq,k,ebpm.fn.l,ebpm.fn.f,res,fix_F,smooth_F,ebps_control){

  # update l
  #l_seq = rowSums(Z)
  l_scale = sum(res$qf$Ef[,k])
  fit = ebpm.fn.l(l_seq,l_scale)
  res$ql$El[,k] = fit$posterior$mean
  res$ql$Elogl[,k] = fit$posterior$mean_log
  res$Hl[k] = calc_H(l_seq,l_scale,fit$log_likelihood,fit$posterior$mean,fit$posterior$mean_log)

  if(!fix_F){
    # update f
    #f_seq = colSums(Z)
    f_scale = sum(res$ql$El[,k])
    if(smooth_F){
      #print(res$gf)
      fit = ebpm.fn.f(f_seq,f_scale,
                      g_init = list(sigma2 = res$gf$sigma2[k]),
                      q_init = list(smooth = res$qf$Elogf_smooth[,k]),
                      general_control = list(maxiter=ebps_control$maxiter,
                                             maxiter_vga = ebps_control$maxiter_vga,
                                             make_power_of_2=ebps_control$make_power_of_2,
                                             vga_tol=ebps_control$vga_tol,
                                             tol = ebps_control$tol),
                      smooth_control = list(wave_trans=ebps_control$wave_trans,
                                            ndwt_method = ebps_control$ndwt_method,
                                            filter.number = ebps_control$filter.number,
                                            family = ebps_control$family,
                                            ebnm_params=ebps_control$ebnm_params,
                                            warmstart=ebps_control$warmstart))
      res$qf$Ef[,k] = fit$posterior$mean
      res$qf$Elogf[,k] = fit$posterior$mean_log
      res$gf$sigma2[k] = fit$fitted_g$sigma2
      res$qf$Ef_smooth[,k] = fit$posterior$mean_smooth
      res$qf$Elogf_smooth[,k] = fit$posterior$mean_log_smooth
    }else{
      fit = ebpm.fn.f(f_seq,f_scale)
      res$qf$Ef[,k] = fit$posterior$mean
      res$qf$Elogf[,k] = fit$posterior$mean_log
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
#'@export
ebpmf_identity_smooth_control_default = function(){
  list(wave_trans='ndwt',
       ndwt_method = "ti.thresh",
       filter.number = 1,
       family = 'DaubExPhase',
       ebnm_params=list(),
       maxiter=1,
       maxiter_vga = 10,
       make_power_of_2='extend',
       vga_tol=1e-3,
       tol = 1e-2,
       warmstart=TRUE)
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

calc_qz = function(n,p,K,ql,qf){
  qz = array(dim = c(n,p,K))
  for(k in 1:K){
    qz[,,k] = outer(ql$Elogl[,k], qf$Elogf[,k], "+")
  }
  return(softmax3d(qz))
}


