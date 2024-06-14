#'@title Fit Sparse Poisson matrix factorization via variational empirical Bayes
#'@param Y count data matrix
#'@param S The known scaling factor matrix, background frequency, default = 11'.
#'@param sigma2 the variance term
#'@param est_sigma2 whether estimate the variance term or fix it
#'@param ebnm.fn see `?flash`.
#'@return fitted object
#'@import flashier
#'@import magrittr
#'@importFrom parallel mclapply
#'@importFrom vebpm pois_mean_GG
#'@export
ebpmf = function(Y,S=NULL,
                  sigma2=NULL,est_sigma2 = TRUE,
                  ebnm.fn = ebnm::ebnm_point_normal,
                  loadings_sign = 0,
                  factors_sign = 0,
                  Kmax=50,
                  var_type='by_col',
                  M_init = NULL,
                  maxiter=100,
                  tol=1e-5,
                  maxiter_backfitting = 1,
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

  if(is.null(S)){
    S = 1
  }
  if(length(S)==1){
    S = matrix(S,nrow=n,ncol=p)
  }

  init_val = ebpmf_init(Y,S,sigma2,
                        var_type,
                        M_init,
                        verbose,
                        n_cores,
                        printevery)
  sigma2 = init_val$sigma2_init
  M = init_val$M_init
  rm(init_val)

  # calc constant for elbo
  const = sum(Y*log(S)) - sum(lfactorial(Y))
  # this is for n=p.
  if(var_type=='by_row'){
    S.dim = 1
  }else if(var_type=='by_col'){
    S.dim = 2
  }else if(var_type=='constant'){
    S.dim = NULL
  }else{
    stop('Non-supported variance type')
  }

  if(verbose){
    cat('running initial flash greedy + backfitting')
    cat('\n')
  }
  #print(sigma2)
  ## fit flashier on M, sigma2 with backfitting
  init.fn = function(f){init.fn.default(f, dim.signs = c(loadings_sign, factors_sign))}
  fit_flash = suppressWarnings(flash.init(M, S = sqrt(sigma2), var.type = NULL, S.dim = S.dim)%>%
    flash.add.greedy(Kmax = Kmax,verbose = verbose_flash,ebnm.fn=ebnm.fn,init.fn=init.fn) %>%
    flash.backfit(verbose = verbose_flash,maxiter=maxiter_backfitting) %>%
    flash.nullcheck(verbose = verbose_flash))

  if(fit_flash$n.factors==0){
    stop('n_factor = 0. How to deal with this issue?')
    }
  K_trace = c(fit_flash$n.factors)
  obj = -Inf
  if(verbose){
    print('Running iterations...')
  }

  for(iter in 1:maxiter){
    # t0 = Sys.time()
    res = vga_pois_solver_mat(M,Y,S,fitted(fit_flash),adjust_var_shape(sigma2,var_type,n,p))
    M = res$M
    V = res$V
    # t1 = Sys.time()
    # run_time_vga[iter] = difftime(t1,t0,units='secs')
    # print(paste('After vga,elbo is',calc_split_PMF_obj_flashier(Y,S,sigma2,M,V,fit_flash,KL_LF,const,var_type)))

    # update sigma2
    if(est_sigma2){
      if(var_type=='constant'){
        sigma2 = (sum(V) +sum(fit_flash$flash.fit$R2))/(n*p)
      }else if(var_type=='by_row'){
        sigma2 = (rowSums(V)+fit_flash$flash.fit$R2)/p
      }else if(var_type=='by_col'){
        sigma2 = (colSums(V)+fit_flash$flash.fit$R2)/n
      }else{
        stop('Non-supported var type')
      }
    }
    #print(paste('After sigma2,elbo is',calc_split_PMF_obj_flashier(Y,S,sigma2,M,V,fit_flash,KL_LF,const,var_type)))
    ## solve flash


    fit_flash = suppressWarnings(flash.init(M, S = sqrt(sigma2), var.type = NULL, S.dim = S.dim)%>%
      flash.init.factors(init = fit_flash,ebnm.fn=ebnm.fn) %>%
      flash.add.greedy(Kmax = Kmax,verbose = verbose_flash,ebnm.fn=ebnm.fn,init.fn=init.fn) %>%
      flash.backfit(verbose = verbose_flash,maxiter = maxiter_backfitting) %>%
      flash.nullcheck(verbose = verbose_flash))


    K_trace[iter+1] = fit_flash$n.factors
    KL_LF = sum(ff.KL(fit_flash$flash.fit,1)) + sum(ff.KL(fit_flash$flash.fit,2))
    #print(paste('After flash,elbo is',calc_split_PMF_obj_flashier(Y,S,sigma2,M,V,fit_flash,KL_LF,const,var_type)))
    # check convergence
    obj[iter + 1] = calc_split_PMF_obj_flashier(Y,S,sigma2,M,V,fit_flash,KL_LF,const,var_type)
    if((obj[iter+1] - obj[iter])/num_points < tol){
      if((obj[iter+1] - obj[iter])<0){
        warning('An iteration decreases ELBO')
      }
      break
    }

    if(verbose){
      if(iter%%printevery==0){
        print(paste('At iter ',iter, ', ELBO=',round(obj[iter+1],log10(1/tol)),sep = ''))
      }
    }

    # save fitted values because the total running time could be very long

    if(iter%%save_fit_every==0){
      saveRDS(list(fit_flash=fit_flash,
                   elbo=obj[length(obj)],
                   K_trace=K_trace,
                   elbo_trace=obj,
                   sigma2 = sigma2,
                   run_time = difftime(Sys.time(),start_time,units='auto')
                   #M=M,V=V,
                   #init_val=init_val,
                   ),
              file=paste(save_fit_path,save_fit_name,'_fit_iter',iter,'.rds',sep=''))
    }


  }

  end_time = Sys.time()
  return(list(fit_flash=fit_flash,
              elbo=obj[length(obj)],
              K_trace=K_trace,
              elbo_trace=obj,
              sigma2 = sigma2,
              run_time = difftime(end_time,start_time,units='auto')
              #M=M,V=V,
              #init_val=init_val
              ))
}

ebpmf_init = function(Y,S,sigma2,
                      var_type,
                      M_init,
                      verbose,
                      n_cores,
                      printevery){
  # we need a. M; b. sigma2
  # case 1: both are given in input.
  if(!is.null(sigma2) & !is.null(M_init)){
    return(list(sigma2_init=sigma2_init,M_init=M_init))
  }

  # case 2: sigma2 is given as input, but M_init is null
  # then we fit poisGG with mean = 0, var =sigma2
  if(is.null(M_init) & !is.null(sigma2)){
    if(verbose){
      cat('Initializing M...')
    }
    # pre-estimate sigma2, assuming LF = 0?.
    if(var_type=='constant'){
      if(verbose){
        cat('Solving VGA...')
      }
      init_val = pois_mean_GG(as.vector(Y),as.vector(S),prior_mean = 0,prior_var = sigma2,tol=1e-3)
      M0 = matrix(init_val$posterior$mean_log,nrow=n,ncol=p)
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
        fit = pois_mean_GG(Y[i,],S[i,],prior_mean = 0,prior_var = sigma2,tol=1e-3)
        return(fit$posterior$mean_log)
      },mc.cores = n_cores)
      M0 = do.call(rbind,init_val)
    }
    if(var_type=='by_col'){
      if(verbose){
        cat('Solving VGA for column 1...')
      }
      init_val = mclapply(1:p,function(i){
        if(verbose){
          if(i%%printevery==0){
            cat(paste(i,'...'))
          }
        }
        fit = pois_mean_GG(Y[,i],S[,i],prior_mean = 0,prior_var = sigma2,tol=1e-3)
        return(fit$posterior$mean_log)
      },mc.cores = n_cores)
      M0 = do.call(cbind,init_val)
    }
    return(list(sigma2_init=sigma2,M_init=M0))
  }
  # case 3: M_init is given, sigma2 is unknown, this case is not interesting and very unrealistic even for development purpose.

  # case 4: both M_init and sigma2 are not given
  if(is.null(sigma2) & is.null(M_init)){
    if(verbose){
      cat('Initializing...')
    }
    # pre-estimate sigma2, assuming LF = 0?.
    if(var_type=='constant'){
      if(verbose){
        cat('Solving VGA...')
      }
      init_val = pois_mean_GG(as.vector(Y),as.vector(S),prior_mean = 0,prior_var = NULL,tol=1e-3)
      sigma2_init = init_val$fitted_g$var
      M0 = matrix(init_val$posterior$mean_log,nrow=n,ncol=p)
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
        fit = pois_mean_GG(Y[i,],S[i,],prior_mean = 0,prior_var = NULL,tol=1e-3)
        return(list(sigma2 = fit$fitted_g$var,mean_log = fit$posterior$mean_log))
      },mc.cores = n_cores)
      sigma2_init = unlist(lapply(init_val,function(fit){fit$sigma2}))
      M0 = do.call(rbind,lapply(init_val,function(fit){fit$mean_log}))
    }
    if(var_type=='by_col'){
      if(verbose){
        cat('Solving VGA for column 1...')
      }
      init_val = mclapply(1:p,function(i){
        if(verbose){
          if(i%%printevery==0){
            cat(paste(i,'...'))
          }
        }
        fit = pois_mean_GG(Y[,i],S[,i],prior_mean = 0,prior_var = NULL,tol=1e-3)
        return(list(sigma2 = fit$fitted_g$var,mean_log = fit$posterior$mean_log))
      },mc.cores = n_cores)
      sigma2_init = unlist(lapply(init_val,function(fit){fit$sigma2}))
      M0 = do.call(cbind,lapply(init_val,function(fit){fit$mean_log}))
    }
    return(list(sigma2_init=sigma2_init,M_init=M0))
  }
}








