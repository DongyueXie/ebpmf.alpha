#'@title Fit empirical Bayes matrix factorization
#'@param Y count data matrix, in sparse format
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
#'@param a0,b0 Inverse-Gamma(a0,b0) prior on sigma2 for regularization.
#'@param cap_var_mean_ratio only update sigma2 when if var/mean > (1+cap_var_mean_ratio). i.e. when overdispersion is low enough, stop updating sigma2 to boost convergence.
#'@param garbage_collection_every perform gc() to reduce memory usage.
#'@return fitted object
#'@import flashier
#'@import magrittr
#'@importFrom parallel mclapply
#'@importFrom vebpm pois_mean_GG
#'@importFrom Matrix Diagonal
#'@importFrom matrixStats colMaxs
#'@importFrom matrixStats rowMaxs
#'@export
ebpmf_log = function(Y,l0=NULL,f0=NULL,
                      var_type='by_col',
                      sigma2=NULL,
                      est_sigma2 = TRUE,
                     M_init = NULL,

                     ebnm.fn = ebnm::ebnm_point_normal,
                     loadings_sign = 0,
                     factors_sign = 0,
                     Kmax_init=30,
                    add_greedy_Kmax = 1,
                    add_greedy_warmstart = TRUE,
                    add_greedy_extrapolate = FALSE,
                    add_greedy_init = 'new_init',
                    add_greedy_every = 1,
                    batch_size = 1e3,
                    maxiter=100,
                    maxiter_backfitting = 1,
                    maxiter_vga = 3,
                    conv_tol=1e-5,
                    init_tol = 1e-5,
                    vga_tol = 1e-3,
                    verbose_flash=0,
                    printevery=10,
                    verbose=FALSE,
                    n_cores = 1,
                    a0 = 1,
                    b0 = 1,
                    cap_var_mean_ratio = 0.1,
                    save_init_val = FALSE,
                    return_sigma2_trace = FALSE,
                    garbage_collection_every = 10,
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

  init_val = ebpmf_log_init(Y,l0,f0,sigma2,
                                var_type,
                                M_init,
                                verbose,
                                n_cores,
                                init_tol,
                                printevery)
  run_time_vga_init = difftime(Sys.time(),start_time,units = 'secs')

  sigma2 = init_val$sigma2_init
  M = init_val$M_init
  if(!save_init_val){
    init_val = NULL
  }
  gc()

  if(is(Y,'sparseMatrix')){
    const = sum(Diagonal(n,log(l0))%*%Y) + sum(Y%*%Diagonal(p,log(f0)))- sum_lfactorial_sparseMat(Y)
  }else{
    const = sum(Y*log(tcrossprod(l0,f0))) - sum(lfactorial(Y))
  }

  if(var_type=='by_row'){
    S.dim = 1
    var_offset_for_obj = p
    var.type = 1
  }else if(var_type=='by_col'){
    S.dim = 2
    var_offset_for_obj = n
    var.type = 2
  }else if(var_type=='constant'){
    S.dim = NULL
    var_offset_for_obj = p*n
    var.type = 0
  }else{
    stop('Non-supported variance type')
  }

  ## split data
  n_batch = ceiling(n/batch_size)
  if(n_batch>1){
    # if(n_batch==1){
    #   batches = list(idx = 1:n)
    # }else{
    #   batches = split(1:n, cut(seq_along(1:n),n_batch , labels = FALSE))
    # }
    batches = split(1:n, cut(seq_along(1:n),n_batch , labels = FALSE))
    # transform Y to a list of sub-Y's
    Y = lapply(batches,function(b){Y[b,]})
    l0 = lapply(batches,function(b){l0[b]})
  }else{
    # To speed up calculation when Y is small and dense Y can be fitted into memory
    Y = as.matrix(Y)
  }

  if(verbose){
    cat('\n')
    cat('running initial flash fit')
    cat('\n')
  }
  if(loadings_sign==0&factors_sign==0){
    # this is faster than the default init method in flash
    init.fn.flash = function(f){init.fn.irlba(f)}
  }else{
    init.fn.flash = function(f){init.fn.default(f, dim.signs = c(loadings_sign, factors_sign))}
  }

  t0 = Sys.time()
  if(is.null(sigma2)){
    fit_flash = suppressWarnings(flash.init(M, S = NULL, var.type = var.type)%>%
                                   flash.add.greedy(Kmax = Kmax_init,verbose = verbose_flash,ebnm.fn=ebnm.fn,init.fn=init.fn.flash,extrapolate = add_greedy_extrapolate) %>%
                                   flash.backfit(verbose = verbose_flash,maxiter = maxiter_backfitting) %>%
                                   flash.nullcheck(verbose = verbose_flash))
    sigma2 = fit_flash$residuals.sd^2
  }else{
    fit_flash = suppressWarnings(flash.init(M, S = sqrt(sigma2), var.type = NULL, S.dim = S.dim)%>%
                                   flash.add.greedy(Kmax = Kmax_init,verbose = verbose_flash,ebnm.fn=ebnm.fn,init.fn=init.fn.flash,extrapolate = add_greedy_extrapolate) %>%
                                   flash.backfit(verbose = verbose_flash,maxiter = maxiter_backfitting) %>%
                                   flash.nullcheck(verbose = verbose_flash))
  }
  rm(M)
  gc()
  run_time_flash_init =  difftime(Sys.time(),t0,units = 'secs')
  if(fit_flash$n.factors==0){
    stop('No structure found in initialization. How to deal with this issue? One way is to reduce sigma2? Or let flash estimate sigma2 instead of fixing it?')
  }

  K_trace = fit_flash$n.factors
  sigma2_trace = sigma2
  K_changed = TRUE
  obj = -Inf

  run_time_vga = c()
  run_time_flash_init_factor = c()
  run_time_flash_greedy = c()
  run_time_flash_backfitting = c()
  run_time_flash_nullcheck = c()

  if(verbose){
    cat('Running iterations...')
    cat('\n')
  }

  for(iter in 1:maxiter){

    #sigma2_old = sigma2

    t0 = Sys.time()
    if(n_batch >1){
      v_sum=0
      sym = 0
      ssexp = 0
      slogv = 0
      for(i_b in 1:n_batch){
        b = batches[[i_b]]
        ## this is for speeding up
        y_b = as.matrix(Y[[i_b]])
        # beta_b = tcrossprod(fit_flash$flash.fit$EF[[1]][b,],fit_flash$flash.fit$EF[[2]])
        # sigma2_b = ifelse(var_type=='by_row',sigma2[b],adjust_var_shape(sigma2,var_type,length(b),p))
        res = vga_pois_solver_mat_newton(fit_flash$flash.fit$Y[b,],
                                         y_b,
                                         tcrossprod(l0[[i_b]],f0),
                                         tcrossprod(fit_flash$flash.fit$EF[[1]][b,],fit_flash$flash.fit$EF[[2]]),
                                         my_ifelse(var_type=='by_row',sigma2[b],adjust_var_shape(sigma2,var_type,length(b),p)),
                                         maxiter=maxiter_vga,tol=vga_tol,return_V=TRUE)

        fit_flash$flash.fit$Y[b,] = res$M
        sym = sym + sum(y_b*res$M)
        ssexp = ssexp + sum(tcrossprod(l0[[i_b]],f0)*exp(res$M+res$V/2))
        slogv = slogv + sum(log(res$V)/2+0.9189385)
        if(var_type=='by_row'){
          v_sum = c(v_sum,rowSums(res$V))
        }else if(var_type=='by_col'){
          v_sum=v_sum+colSums(res$V)
        }else if(var_type=='constant'){
          v_sum=v_sum+sum(res$V)
        }
      }
      rm(y_b)
      #rm(beta_b)
      #rm(sigma2_b)
      rm(res)
      if(var_type=='by_row'){
        v_sum = v_sum[-1]
      }

      if(est_sigma2){
        sigma2 = ebpmf_update_sigma2(fit_flash,sigma2,v_sum,var_type,cap_var_mean_ratio,a0,b0,n,p)
      }

    }else{
      res = vga_pois_solver_mat_newton(fit_flash$flash.fit$Y,Y,tcrossprod(l0,f0),fitted(fit_flash),
                                       adjust_var_shape(sigma2,var_type,n,p),
                                       maxiter = maxiter_vga,
                                       tol=vga_tol,return_V = TRUE)
      fit_flash$flash.fit$Y = res$M
      sym = sum(Y*res$M)
      ssexp = sum(tcrossprod(l0,f0)*exp(res$M+res$V/2))
      slogv = sum(log(res$V)/2+0.9189385)
      if(var_type=='constant'){
        v_sum =sum(res$V)
      }else if(var_type=='by_col'){
        v_sum =colSums(res$V)
      }else if(var_type=='by_row'){
        v_sum =rowSums(res$V)
      }


      if(est_sigma2){
        sigma2=ebpmf_update_sigma2(fit_flash,sigma2,v_sum,var_type,cap_var_mean_ratio,a0,b0,n,p)
      }
      rm(res)
    }

    t1 = Sys.time()
    run_time_vga[iter] = difftime(t1,t0,units='secs')
    if(return_sigma2_trace){
      sigma2_trace = rbind(sigma2_trace,sigma2)
    }

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

    # check convergence
    obj[iter + 1] = calc_ebpmf_log_obj(n,p,sym,ssexp,slogv,v_sum,sigma2,fit_flash$flash.fit$R2,KL_LF,const,var_offset_for_obj,a0,b0)
    if((obj[iter + 1]-obj[iter])/num_points< conv_tol){
      break
    }
    if(verbose){
      if(iter%%printevery==0){
        cat(paste('iter ',iter, ', elbo=',round(obj[iter+1],log10(1/conv_tol)),", K=",fit_flash$n.factors,sep = ''))
        cat('\n')
      }
    }


    if(iter%%save_fit_every==0){
      saveRDS(list(fit_flash=list(L.pm=fit_flash$L.pm,F.pm = fit_flash$F.pm,pve = fit_flash$pve),
                   elbo=obj[length(obj)],
                   K_trace=K_trace,
                   elbo_trace=obj,
                   sigma2 = sigma2,
                   sigma2_trace = sigma2_trace,
                   run_time = difftime(Sys.time(),start_time,units='auto'),
                   run_time_break_down = list(run_time_vga_init = run_time_vga_init,
                                              run_time_flash_init = run_time_flash_init,
                                              run_time_vga = run_time_vga,
                                              run_time_flash_init_factor = run_time_flash_init_factor,
                                              run_time_flash_greedy = run_time_flash_greedy,
                                              run_time_flash_backfitting = run_time_flash_backfitting,
                                              run_time_flash_nullcheck = run_time_flash_nullcheck)),
              file=paste(save_fit_path,save_fit_name,'_iter',iter,'.rds',sep=''))
    }

    if(iter%%garbage_collection_every==0){
      gc()
    }


  }

  end_time = Sys.time()
  return(list(fit_flash=fit_flash,
              elbo=obj[length(obj)],
              K_trace=K_trace,
              elbo_trace=obj,
              sigma2 = sigma2,
              sigma2_trace = sigma2_trace,
              init_val=init_val,
              run_time = difftime(end_time,start_time,units='auto'),
              run_time_break_down = list(run_time_vga_init = run_time_vga_init,
                                         run_time_flash_init = run_time_flash_init,
                                         run_time_vga = run_time_vga,
                                         run_time_flash_init_factor = run_time_flash_init_factor,
                                         run_time_flash_greedy = run_time_flash_greedy,
                                         run_time_flash_backfitting = run_time_flash_backfitting,
                                         run_time_flash_nullcheck = run_time_flash_nullcheck)))
}

#'@title Calc elbo
calc_ebpmf_log_obj = function(n,p,sym,ssexp,slogv,sv,sigma2,R2,KL_LF,const,ss,a0,b0){
  val = sym - ssexp + slogv + 0.5*n*p - sum(ss*log(2*pi*sigma2)/2)- sum(sv/2/sigma2) - sum(R2/2/sigma2) + const+ KL_LF - sum((a0+1)*log(sigma2)) - sum(b0/sigma2)
  # val = sum(Y*M - S*exp(M+V/2)   + log(2*pi*V)/2 + 0.5 ) - sum(ss*log(2*pi*sigma2)/2)- sum(sv/2/sigma2) - sum(R2/2/sigma2) + const+ KL_LF
  return(val)
}

#'@title update sigma2
ebpmf_update_sigma2 = function(fit_flash,sigma2,v_sum,var_type,cap_var_mean_ratio,a0,b0,n,p){
  if(var_type=='constant'){
    if(cap_var_mean_ratio>0){
      if(((exp(sigma2)-1)*exp(max(fitted(fit_flash))))>cap_var_mean_ratio){
        sigma2 = ((v_sum +sum(fit_flash$flash.fit$R2))/2+b0)/(n*p/2+a0+1)
      }
    }else{
      sigma2 = ((v_sum +sum(fit_flash$flash.fit$R2))/2+b0)/(n*p/2+a0+1)
    }
  }else if(var_type=='by_row'){
    if(cap_var_mean_ratio>0){
      update_idx = which(((exp(sigma2)-1)*exp(rowMaxs(fitted(fit_flash))))>cap_var_mean_ratio)
      if(!is.null(update_idx)){
        sigma2[update_idx] = (((v_sum+fit_flash$flash.fit$R2)/2+b0)/(p/2+a0+1))[update_idx]
      }
    }else{
      sigma2 = ((v_sum+fit_flash$flash.fit$R2)/2+b0)/(p/2+a0+1)
    }
  }else if(var_type=='by_col'){
    if(cap_var_mean_ratio>0){
      update_idx = which(((exp(sigma2)-1)*exp(colMaxs(fitted(fit_flash))))>cap_var_mean_ratio)
      if(!is.null(update_idx)){
        sigma2[update_idx] = (((v_sum+fit_flash$flash.fit$R2)/2+b0)/(n/2+a0+1))[update_idx]
      }
    }else{
      sigma2 = ((v_sum+fit_flash$flash.fit$R2)/2+b0)/(n/2+a0+1)
    }
  }else{
    stop('Non-supported var type')
  }
  return(sigma2)
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
    return(rowSums(Sigma2/(Sigma2*Y-M+Beta+1)))
  }
}

#'@title calculate sum of log factorial of all elements of a sparse Matrix
sum_lfactorial_sparseMat = function(Y){
  sum(lfactorial(Y@x))
}

#'@title adjust var shape for vga. For computation purpose.
adjust_var_shape = function(sigma2,var_type,n,p){
  if(var_type=='constant'){
    sigma2 = matrix(sigma2,nrow=n,ncol=p)
  }else if(var_type=='by_row'){
    sigma2 = matrix(sigma2,nrow=n,ncol=p,byrow = F)
  }else if(var_type=='by_col'){
    sigma2 = matrix(sigma2,nrow=n,ncol=p,byrow = T)
  }else{
    stop('Non-supported var type')
  }
  sigma2
}

my_ifelse = function(test,yes,no){
  if(test){
    return(tes)
  }else{
    return(no)
  }
}

