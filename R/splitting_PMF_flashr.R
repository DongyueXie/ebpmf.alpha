#'@title Fit Sparse Poisson matrix factorization
#'@param Y count data matrix
#'@param S The known scaling factor matrix, background frequency.
#'@param sigma2 the variance term
#'@param est_sigma2 whether estimate the variance term or fix it
#'@return fitted object

splitting_PMF_flashr = function(Y,S,sigma2=NULL,est_sigma2 = TRUE,
                         Kmax=10,var_type='by_row',
                         M_init = NULL,
                         maxiter=100,tol=0.01,
                         verbose_flash=FALSE,
                         printevery=10){

  n = nrow(Y)
  p = ncol(Y)

  if(is.null(S)){
    S = 1
  }
  if(length(S)==1){
    S = matrix(S,nrow=n,ncol=p)
  }

  if(is.null(sigma2)|is.null(M_init)){
    # pre-estimate sigma2, assuming LF = 0?.
    if(var_type=='constant'){
      init_val = pois_mean_GG(as.vector(Y),as.vector(S),prior_mean = 0,prior_var = NULL,tol=1e-3)
      sigma2_init = init_val$fitted_g$var
      M0 = matrix(init_val$posterior$mean_log,nrow=n,ncol=p)
    }
    if(var_type=='by_row'){
      init_val = lapply(1:n,function(i){
        return(pois_mean_GG(Y[i,],S[i,],prior_mean = 0,prior_var = NULL,tol=1e-3))
      })
      sigma2_init = unlist(lapply(init_val,function(fit){fit$fitted_g$var}))
      M0 = do.call(rbind,lapply(init_val,function(fit){fit$posterior$mean_log}))
    }
    if(var_type=='by_col'){
      init_val = lapply(1:p,function(i){
        return(pois_mean_GG(Y[,i],S[,i],prior_mean = 0,prior_var = NULL,tol=1e-3))
      })
      sigma2_init = unlist(lapply(init_val,function(fit){fit$fitted_g$var}))
      M0 = do.call(cbind,lapply(init_val,function(fit){fit$posterior$mean_log}))
    }
  }
  if(is.null(sigma2)){
    sigma2 = sigma2_init
    est_sigma2 = TRUE
    rm(sigma2_init)
  }
  if(is.null(M_init)){
    M = M0
    rm(M0)
    rm(init_val)
  }else{
    M = M_init
  }

  const = sum(Y*log(S)) - sum(lfactorial(Y))

  sigma2 = adjust_var_shape(sigma2,var_type,n,p)
  fit_flash = flash(flash_set_data(M,S=sqrt(sigma2)), Kmax=Kmax,var_type = 'zero',verbose =verbose_flash)

  KL_LF = sum(unlist(fit_flash$fit$KL_f)) + sum(unlist(fit_flash$fit$KL_l))
  V = matrix(1/n,nrow=n,ncol=p)
  obj = calc_split_PMF_obj_flashr(Y,S,sigma2,M,V,fit_flash,KL_LF,const)

  start_time = Sys.time()

  for(iter in 1:maxiter){

    res = vga_pois_solver_mat(M,Y,S,fit_flash$fitted_values,sigma2)
    M = res$M
    V = res$V
    print(paste('After vga,elbo is',calc_split_PMF_obj_flashr(Y,S,sigma2,M,V,fit_flash,KL_LF,const)))

    fit_flash = flash(flash_set_data(M,S=sqrt(sigma2)), Kmax=Kmax,var_type = 'zero',verbose =verbose_flash)
    KL_LF = sum(unlist(fit_flash$fit$KL_f)) + sum(unlist(fit_flash$fit$KL_l))
    print(paste('After flash,elbo is',calc_split_PMF_obj_flashr(Y,S,sigma2,M,V,fit_flash,KL_LF,const)))


    # update sigma2
    if(est_sigma2){
      #R2 = flashr:::flash_get_R2(flash_set_data(M),fit_flash$fit)
      R2 = calc_ER2(M,fit_flash)
      if(var_type=='constant'){
        sigma2 = mean(R2+V)
        sigma2 = matrix(sigma2,nrow=n,ncol=p)
      }else if(var_type=='by_row'){
        sigma2 = rowMeans(V+R2)
        sigma2 = matrix(sigma2,nrow=n,ncol=p,byrow = F)
      }else if(var_type=='by_col'){
        sigma2 = colMeans(V+R2)
        sigma2 = matrix(sigma2,nrow=n,ncol=p,byrow = T)
      }else{
        stop('Non-supported var type')
      }
    }

    print(paste('After sigma2,elbo is',calc_split_PMF_obj_flashr(Y,S,sigma2,M,V,fit_flash,KL_LF,const)))

    # check convergence
    obj[iter + 1] = calc_split_PMF_obj_flashr(Y,S,sigma2,M,V,fit_flash,KL_LF,const)
    if((obj[iter+1] - obj[iter]) < tol){
      if((obj[iter+1] - obj[iter])<0){
        warning('An iteration decreases ELBO')
      }
      break
    }

    if(iter%%printevery==0){
      print(paste('At iter',iter, 'ELBO=',obj[iter+1]))
    }

  }

  end_time = Sys.time()
  return(list(fit_flash=fit_flash,elbo=obj[length(obj)],eblo_trace=obj,
              sigma2 = sigma2,run_time = difftime(end_time,start_time,units='auto'),M=M,V=V))
}



calc_split_PMF_obj_flashr = function(Y,S,sigma2,M,V,fit_flash,KL_LF,const){
  sum(Y*M - S*exp(M+V/2) - log(2*pi*sigma2)/2 - (V+calc_ER2(M,fit_flash))/2/sigma2 + log(2*pi*V)/2 + 0.5 ) + const+ KL_LF
}

calc_ER2 = function(Y,f){
  if (is.null(f$fit$EL)) {
    return(Y^2)
  }
  else {
    return((Y - f$fitted_values)^2 + tcrossprod(f$fit$EL2, f$fit$EF2) - tcrossprod(f$fit$EL^2,f$fit$EF^2))
  }
}

