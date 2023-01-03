#'@title Fit Sparse Poisson matrix factorization low memory version
#'@param Y count data matrix
#'@param l0,f0 The known scaling loadings and factors.
#'@param sigma2 the variance term
#'@param est_sigma2 whether estimate the variance term or fix it at sigma2
#'@param ebnm.fn see `?flash`.
#'@param loadings_sign,factors_sign see `?init.fn.default`, must match ebnm.fn
#'@param Kmax_init the Kmax in the first flash fit
#'@param add_greedy_Kmax THe Kmax in add_greedy in iterations
#'@param add_greedy_warmstart,add_greedy_extrapolate see `?flash.add.greedy`
#'@param add_greedy_init either 'previous_init' or "new_init"
#'@param add_greedy_every perform flash add greedy every `add_greedy_every` iterations.
#'@param maxiter,maxiter_backfitting,maxiter_vga max iterations for the splitting, backfitting, vga.
#'@param conv_tol,init_tol,vga_tol tolerance for convergence, initialization vga fit, and vga fit in iterations
#'@param batch_size reduce memory usage for vga step by looping subsets of dataset.
#'@return fitted object
#'@import flashier
#'@import magrittr
#'@importFrom parallel mclapply
#'@importFrom vebpm pois_mean_GG
#'@importFrom Matrix Diagonal
#'@export
splitting_PMF_flashier_low_memory = function(Y,l0=NULL,f0=NULL,
                                  var_type='by_col',
                                  sigma2=NULL,
                                  est_sigma2 = TRUE,
                                  ebnm.fn = ebnm::ebnm_point_normal,
                                  loadings_sign = 0,
                                  factors_sign = 0,
                                  Kmax_init=50,
                                  add_greedy_Kmax = 1,
                                  add_greedy_warmstart = TRUE,
                                  add_greedy_extrapolate = FALSE,
                                  add_greedy_init = 'previous_init',
                                  add_greedy_every = 1,
                                  batch_size = 1e3,
                                  M_init = NULL,
                                  maxiter=100,
                                  maxiter_backfitting = 1,
                                  maxiter_vga = 1,
                                  conv_tol=1e-5,
                                  init_tol = 1e-5,
                                  vga_tol = 1e-5,
                                  verbose_flash=0,
                                  printevery=10,
                                  verbose=FALSE,
                                  n_cores = 1,
                                  save_fit_every = Inf,
                                  save_fit_path = NULL,
                                  save_fit_name = NULL){

  start_time = Sys.time()

  n = nrow(Y)
  p = ncol(Y)
  num_points = n*p

  if(is.null(l0)){
    l0 = 1
  }
  if(is.null(f0)){
    f0 = 1
  }
  if(length(l0)==1){
    l0 = rep(l0,n)
  }
  if(length(f0)==1){
    f0 = rep(f0,p)
  }

  if(is.null(sigma2)|is.null(M_init)){
    if(verbose){
      cat('Initializing...')
    }
    ## pre-estimate sigma2, assuming LF = 0?.
    if(var_type=='constant'){
      if(verbose){
        cat('Solving VGA constant...')
      }
      init_val = suppressWarnings(pois_mean_GG(as.vector(Y),as.vector(tcrossprod(l0,f0)),prior_mean = 0,prior_var = NULL,tol=init_tol))
      sigma2_init = init_val$fitted_g$var
      M = matrix(init_val$posterior$mean_log,nrow=n,ncol=p)
      rm(init_val)
      # V = matrix(init_val$posterior$var_log,nrow=n,ncol=p)
    }
    if(var_type=='by_row'){
      if(verbose){
        cat('Solving VGA for row 1...')
      }
      init_val = mclapply(1:n,function(i){
        if(verbose){
          if(i%%printevery==0){
            cat(paste(i,'...'))
          }
        }
        fit = suppressWarnings(pois_mean_GG(Y[i,],l0[i]*f0,prior_mean = 0,prior_var = NULL,tol=init_tol))
        return(list(sigma2 = fit$fitted_g$var,mean_log = fit$posterior$mean_log))
      },mc.cores = n_cores)
      sigma2_init = unlist(lapply(init_val,function(fit){fit$sigma2}))
      M = do.call(rbind,lapply(init_val,function(fit){fit$mean_log}))
      rm(init_val)
    }
    if(var_type=='by_col'){
      if(verbose){
        cat('Solving VGA for column 1...')
      }
      if(n_cores>1){
        init_val = mclapply(1:p,function(i){
          if(verbose){
            if(i%%printevery==0){
              cat(paste(i,'...'))
            }
          }
          fit = suppressWarnings(pois_mean_GG(Y[,i],l0*f0[i],prior_mean = 0,prior_var = NULL,tol=init_tol))
          return(list(sigma2 = fit$fitted_g$var,mean_log = fit$posterior$mean_log))
        },mc.cores = n_cores)
        sigma2_init = unlist(lapply(init_val,function(fit){fit$sigma2}))
        M = do.call(cbind,lapply(init_val,function(fit){fit$mean_log}))
        rm(init_val)
      }else{
        M = matrix(nrow=n,ncol=p)
        sigma2_init = rep(0,p)
        for(j in 1:p){
          if(verbose){
            if(j%%printevery==0){
              cat(paste(j,'...'))
            }
          }
          fit = suppressWarnings(pois_mean_GG(Y[,j],l0*f0[j],prior_mean = 0,prior_var = NULL,tol=init_tol))
          M[,j] = fit$posterior$mean_log
          sigma2_init[j] = fit$fitted_g$var
          rm(fit)
        }
      }
    }
    gc()
  }
  run_time_vga_init = difftime(Sys.time(),start_time,units = 'secs')

  if(is.null(sigma2)){
    sigma2 = sigma2_init
    est_sigma2 = TRUE
    rm(sigma2_init)
  }
  if(!is.null(M_init)){
    M = M_init
    rm(M_init)
  }

   const = sum(Diagonal(n,log(l0))%*%Y) + sum(Y%*%Diagonal(p,log(f0)))- sum_lfactorial_sparseMat(Y)

  if(var_type=='by_row'){
    S.dim = 1
  }else if(var_type=='by_col'){
    S.dim = 2
  }else if(var_type=='constant'){
    S.dim = NULL
  }else{
    stop('Non-supported variance type')
  }

  ## split data
  n_batch = ceiling(n/batch_size)
  if(n_batch>1){
    batches = split(1:n, cut(seq_along(1:n),n_batch , labels = FALSE))
    # transform Y to a list of sub-Y's
    Y = lapply(batches,function(b){Y[b,]})
    l0 = lapply(batches,function(b){l0[b]})
  }else{
    # To speed up calculation when Y is small and dense Y can be fitted into memory
    Y = as.matrix(Y)
  }

  if(verbose){
    cat('running initial flash fit')
    cat('\n')
  }
  if(loadings_sign==0&factors_sign==0){
    init.fn.flash = function(f){init.fn.irlba(f)}
  }else{
    init.fn.flash = function(f){init.fn.default(f, dim.signs = c(loadings_sign, factors_sign))}
  }

  t0 = Sys.time()
  fit_flash = suppressWarnings(flash.init(M, S = sqrt(sigma2), var.type = NULL, S.dim = S.dim)%>%
                                 flash.add.greedy(Kmax = Kmax_init,verbose = verbose_flash,ebnm.fn=ebnm.fn,init.fn=init.fn.flash,extrapolate = add_greedy_extrapolate) %>%
                                 flash.backfit(verbose = verbose_flash,maxiter = maxiter_backfitting) %>%
                                 flash.nullcheck(verbose = verbose_flash))
  rm(M)
  gc()
  run_time_flash_init =  difftime(Sys.time(),t0,units = 'secs')
  if(fit_flash$n.factors==0){
    stop('No structure found in initialization. How to deal with this issue?')
  }

  K_trace = fit_flash$n.factors
  K_changed = TRUE
  obj = -Inf

  run_time_vga = c()
  run_time_flash_init_factor = c()
  run_time_flash_greedy = c()
  run_time_flash_backfitting = c()
  run_time_flash_nullcheck = c()

  if(verbose){
    print('Running iterations...')
  }

  # KL_LF = sum(ff.KL(fit_flash$flash.fit,1)) + sum(ff.KL(fit_flash$flash.fit,2))

  for(iter in 1:maxiter){

    sigma2_old = sigma2

    t0 = Sys.time()
    if(n_batch > 1){
      v_sum=0
      sym = 0
      ssexp = 0
      slogv = 0
      for(i_b in 1:n_batch){
        b = batches[[i_b]]
        # this is for speeding up
        y_b = as.matrix(Y[[i_b]])
        beta_b = tcrossprod(fit_flash$flash.fit$EF[[1]][b,],fit_flash$flash.fit$EF[[2]])
        sigma2_b = ifelse(var_type=='by_row',sigma2[b],adjust_var_shape(sigma2,var_type,length(b),p))
        res = vga_pois_solver_mat_newton(fit_flash$flash.fit$Y[b,],
                                                               y_b,tcrossprod(l0[[i_b]],f0),
                                                                beta_b,
                                                                sigma2_b,
                                                                maxiter=maxiter_vga,tol=vga_tol,return_V=TRUE)
        fit_flash$flash.fit$Y[b,] = res$M
        sym = sym + sum(y_b*res$M)
        ssexp = ssexp + sum(tcrossprod(l0[[i_b]],f0)*exp(res$M+res$V/2))
        slogv = slogv + sum(log(res$V)/2+0.9189385)
        if(est_sigma2){
          if(var_type=='by_row'){
            # sigma2[b] = (V_M(fit_flash$flash.fit$Y[b,],y_b,beta_b,sigma2_b,var_type='by_row')+fit_flash$flash.fit$R2[b])/p
            v_sum = c(v_sum,rowSums(res$V))
          }else if(var_type=='by_col'){
            v_sum=v_sum+colSums(res$V)
          }else if(var_type=='constant'){
            v_sum=v_sum+sum(res$V)
          }
        }
      }
      rm(y_b)
      rm(beta_b)
      rm(sigma2_b)
      rm(res)
      if(est_sigma2){
        if(var_type=='constant'){
          sigma2 = (v_sum+sum(fit_flash$flash.fit$R2))/(n*p)
        }else if(var_type=='by_col'){
          sigma2 = (v_sum+fit_flash$flash.fit$R2)/n
        }else if(var_type=='by_row'){
          v_sum = v_sum[-1]
          sigma2 = (v_sum+fit_flash$flash.fit$R2)/p
        }
      }
    }else{
      # #beta_b = tcrossprod(fit_flash$flash.fit$EF[[1]],fit_flash$flash.fit$EF[[2]])
      # beta_b = fitted(fit_flash)
      # sigma2_b = ifelse(var_type=='by_row',sigma2,adjust_var_shape(sigma2,var_type,n,p))
      # fit_flash$flash.fit$Y = vga_pois_solver_mat_newton(fit_flash$flash.fit$Y,Y,tcrossprod(l0,f0),beta_b,
      #                                                    sigma2_b,
      #                                                    maxiter = maxiter_vga,
      #                                                    tol=vga_tol,
      #                                                    return_V=FALSE)
      # # update sigma2
      # if(est_sigma2){
      #   if(var_type=='constant'){
      #     sigma2 = (V_M(fit_flash$flash.fit$Y,Y,beta_b,sigma2_b,var_type='constant')+sum(fit_flash$flash.fit$R2))/(n*p)
      #   }else if(var_type=='by_row'){
      #     sigma2 = (V_M(fit_flash$flash.fit$Y,Y,beta_b,sigma2_b,var_type='by_row')+fit_flash$flash.fit$R2)/p
      #   }else if(var_type=='by_col'){
      #     sigma2 = (V_M(fit_flash$flash.fit$Y,Y,beta_b,sigma2_b,var_type='by_col')+fit_flash$flash.fit$R2)/n
      #   }else{
      #     stop('Non-supported var type')
      #   }
      # }
      # rm(beta_b)
      # rm(sigma2_b)
    }

    t1 = Sys.time()
    run_time_vga[iter] = difftime(t1,t0,units='secs')



    gc()
    # if(iter>1){
    #   print(paste('iteration:',iter))
    #   print(paste('vga,elbo is',round(calc_split_PMF_obj_flashier(Y,S,sigma2,M,V,fit_flash,KL_LF,const,var_type),3)))
    # }
    # if(iter>1){
    #   print(paste('sigma2,elbo is',round(calc_split_PMF_obj_flashier(Y,S,sigma2,M,V,fit_flash,KL_LF,const,var_type),3)))
    # }

    ## solve flash
    ## To timing the operations, I separate flash fits:
    t0 = Sys.time()
    fit_flash = flash.init(fit_flash$flash.fit$Y, S = sqrt(sigma2), var.type = NULL, S.dim = S.dim) %>%
      flash.init.factors(init = fit_flash,ebnm.fn=ebnm.fn)
    t2 = Sys.time()
    run_time_flash_init_factor[iter] = difftime(t2,t0,units='secs')

    if(iter%%add_greedy_every==0 & fit_flash$n.factors < Kmax_init){
      if(add_greedy_init=='previous_init'){
        if(K_changed){
          init_vals = do.call(init.fn.flash,list(fit_flash$flash.fit))
        }
        fit_flash$flash.fit$init_vals = init_vals
        fit_flash = flash.add.greedy(fit_flash, Kmax = 1,verbose = verbose_flash,
                                     ebnm.fn=ebnm.fn,init.fn = init.fn.fix,
                                     warmstart = add_greedy_warmstart,
                                     extrapolate = add_greedy_extrapolate)
      }else if(add_greedy_init=='new_init'){
        fit_flash = flash.add.greedy(fit_flash, Kmax = add_greedy_Kmax,verbose = verbose_flash,
                                     ebnm.fn=ebnm.fn,init.fn = init.fn.flash,
                                     warmstart = add_greedy_warmstart,
                                     extrapolate = add_greedy_extrapolate)
      }
      K_changed = (fit_flash$n.factors != K_trace[iter])
    }
    t3 = Sys.time()
    run_time_flash_greedy[iter] = difftime(t3,t2,units='secs')

    fit_flash = suppressWarnings(flash.backfit(fit_flash, verbose = verbose_flash,maxiter = maxiter_backfitting))
    t4 = Sys.time()
    run_time_flash_backfitting[iter] = difftime(t4,t3,units='secs')

    fit_flash = flash.nullcheck(fit_flash, verbose = verbose_flash)
    t5 = Sys.time()
    run_time_flash_nullcheck[iter] = difftime(t5,t4,units='secs')

    K_trace[iter+1] = fit_flash$n.factors
    KL_LF = sum(ff.KL(fit_flash$flash.fit,1)) + sum(ff.KL(fit_flash$flash.fit,2))
    # if(iter>1){
    #   print(paste('flash,elbo is',round(calc_split_PMF_obj_flashier(Y,S,sigma2,M,V,fit_flash,KL_LF,const,var_type),3)))
    #   cat("----------------------------")
    #   cat("\n")
    # }

    # check convergence
    obj[iter + 1] = calc_split_PMF_obj_flashier_low_memory(n,p,sym,ssexp,slogv,v_sum,sigma2,fit_flash$flash.fit$R2,KL_LF,const,var_type)
    if((obj[iter + 1]-obj[iter])< conv_tol){
      break
    }
    if(verbose){
      if(iter%%printevery==0){
        print(paste('iter ',iter, ', elbo=',round(obj[iter+1],log10(1/conv_tol)),", K=",fit_flash$n.factors,sep = ''))
      }
    }


    gc()

  }

  end_time = Sys.time()
  return(list(fit_flash=fit_flash,
              elbo=obj[length(obj)],
              K_trace=K_trace,
              elbo_trace=obj,
              sigma2 = sigma2,
              run_time = difftime(end_time,start_time,units='auto'),
              run_time_break_down = list(run_time_vga_init = run_time_vga_init,
                                         run_time_flash_init = run_time_flash_init,
                                         run_time_vga = run_time_vga,
                                         run_time_flash_init_factor = run_time_flash_init_factor,
                                         run_time_flash_greedy = run_time_flash_greedy,
                                         run_time_flash_backfitting = run_time_flash_backfitting,
                                         run_time_flash_nullcheck = run_time_flash_nullcheck)))
}


calc_split_PMF_obj_flashier_low_memory = function(n,p,sym,ssexp,slogv,sv,sigma2,R2,KL_LF,const,var_type){
  # R2 = fit_flash$flash.fit$R2
  # n = nrow(Y)
  # p = ncol(Y)
  if(var_type=='by_row'){
    #sv = rowSums(V)
    ss = p
  }else if(var_type=='by_col'){
    #sv = colSums(V)
    ss = n
  }else if(var_type=='constant'){
    #sv = sum(V)
    ss = n*p
  }else{
    stop('Non-supported var type')
  }
  val = sym - ssexp + slogv + 0.5*n*p - sum(ss*log(2*pi*sigma2)/2)- sum(sv/2/sigma2) - sum(R2/2/sigma2) + const+ KL_LF
  #val = sum(Y*M - S*exp(M+V/2)   + log(2*pi*V)/2 + 0.5 ) - sum(ss*log(2*pi*sigma2)/2)- sum(sv/2/sigma2) - sum(R2/2/sigma2) + const+ KL_LF
  return(val)
}

# V_M = function(M,Y,EL,EF,sigma2,var_type='by_col'){
#   if(var_type=='by_col'){
#     p = ncol(M)
#     return(colSums((1/(Y%*%Diagonal(p,sigma2)+tcrossprod(EL,EF)-M+1))%*%Diagonal(p,sigma2)))
#   }
#   if(var_type=='by_row'|var_type=='constant'){
#     return(rowSums(sigma2/(simga2*Y-M+tcrossprod(EL,EF)+1)))
#   }
# }

#'@title Calc V given M
V_M = function(M,Y,Beta,Sigma2,var_type='by_col'){
  if(var_type=='by_col'){
    return(colSums(Sigma2/(Y*Sigma2+Beta-M+1)))
  }
  if(var_type=='by_row'|var_type=='constant'){
    return(rowSums(sigma2/(simga2*Y-M+Beta+1)))
  }
}

#'@title calculate sum of log factorial of all elements of a sparse Matrix
sum_lfactorial_sparseMat = function(Y){
  sum(lfactorial(Y@x))
}
