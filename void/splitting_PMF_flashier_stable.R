#'@title Fit Sparse Poisson matrix factorization
#'@param Y count data matrix
#'@param S The known scaling factor matrix, background frequency.
#'@param sigma2 the variance term
#'@param est_sigma2 whether estimate the variance term or fix it
#'@param ebnm.fn see `?flash`.
#'@return fitted object
#'@import flashier
#'@import magrittr
#'@importFrom parallel mclapply
#'@importFrom vebpm pois_mean_GG
#'@export
splitting_PMF_flashier = function(Y,S=NULL,
                                  sigma2=NULL,est_sigma2 = TRUE,
                                  ebnm.fn = ebnm::ebnm_point_normal,
                                  loadings_sign = 0,
                                  factors_sign = 0,
                                  Kmax=50,
                                  var_type='by_col',
                                  M_init = NULL,
                                  maxiter=100,
                                  conv_tol=1e-5,
                                  init_tol = 1e-5,
                                  vga_tol = 1e-5,
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

  if(is.null(sigma2)|is.null(M_init)){
    if(verbose){
      cat('Initializing...')
    }
    # pre-estimate sigma2, assuming LF = 0?.
    if(var_type=='constant'){
      if(verbose){
        cat('Solving VGA...')
      }
      init_val = pois_mean_GG(as.vector(Y),as.vector(S),prior_mean = 0,prior_var = NULL,tol=init_tol)
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
        fit = pois_mean_GG(Y[i,],S[i,],prior_mean = 0,prior_var = NULL,tol=init_tol)
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
        fit = pois_mean_GG(Y[,i],S[,i],prior_mean = 0,prior_var = NULL,tol=init_tol)
        return(list(sigma2 = fit$fitted_g$var,mean_log = fit$posterior$mean_log))
      },mc.cores = n_cores)
      sigma2_init = unlist(lapply(init_val,function(fit){fit$sigma2}))
      M0 = do.call(cbind,lapply(init_val,function(fit){fit$mean_log}))
    }
  }
  run_time_vga_init = difftime(Sys.time(),start_time,units = 'secs')
  init_val = list()
  if(is.null(sigma2)){
    sigma2 = sigma2_init
    est_sigma2 = TRUE
    init_val$sigma2_init = sigma2_init
    rm(sigma2_init)
  }
  if(is.null(M_init)){
    M = M0
    init_val$M_init = M0
    rm(M0)
  }else{
    M = M_init
  }

  const = sum(Y*log(S)) - sum(lfactorial(Y))
  ## fit flashier on M, sigma2 with backfitting
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

  init.fn = function(f){init.fn.default(f, dim.signs = c(loadings_sign, factors_sign))}

  #print(sigma2)
  fit_flash = flash.init(M, S = sqrt(sigma2), var.type = NULL, S.dim = S.dim)%>%
    flash.add.greedy(Kmax = Kmax,verbose = verbose_flash,ebnm.fn=ebnm.fn,init.fn=init.fn) %>%
    flash.backfit(verbose = verbose_flash,maxiter = maxiter_backfitting) %>%
    flash.nullcheck(verbose = verbose_flash)


  if(fit_flash$n.factors==0){
    stop('No structure found in initialization. How to deal with this issue?')
  }

  # KL_LF = sum(ff.KL(fit_flash$flash.fit,1)) + sum(ff.KL(fit_flash$flash.fit,2))
  # V = matrix(1/n,nrow=n,ncol=p)
  # obj = calc_split_PMF_obj_flashier(Y,S,sigma2,M,V,fit_flash,KL_LF,const,var_type)

  K_trace = fit_flash$n.factors
  obj = -Inf



  if(verbose){
    print('Running iterations...')
  }

  run_time_vga = c()
  run_time_flash_init = c()
  run_time_flash_init_factor = c()
  run_time_flash_greedy = c()
  run_time_flash_backfitting = c()
  run_time_flash_nullcheck = c()


  for(iter in 1:maxiter){

    t0 = Sys.time()
    res = vga_pois_solver_mat(M,Y,S,fitted(fit_flash),adjust_var_shape(sigma2,var_type,n,p),tol=vga_tol)
    M = res$M
    V = res$V
    t1 = Sys.time()
    run_time_vga[iter] = difftime(t1,t0,units='secs')
    #print(paste('After vga,elbo is',calc_split_PMF_obj_flashier(Y,S,sigma2,M,V,fit_flash,KL_LF,const,var_type)))

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

    #print(M[1:5,1:5])
    ## To timeing the operations, I seperate flash fits:

    t0 = Sys.time()
    fit_flash_new = flash.init(M, S = sqrt(sigma2), var.type = NULL, S.dim = S.dim)
    t1 = Sys.time()
    run_time_flash_init[iter] = difftime(t1,t0,units='secs')

    fit_flash = flash.init.factors(fit_flash_new,init = fit_flash,ebnm.fn=ebnm.fn)
    t2 = Sys.time()
    run_time_flash_init_factor[iter] = difftime(t2,t1,units='secs')

    fit_flash = flash.add.greedy(fit_flash, Kmax = Kmax,verbose = verbose_flash,ebnm.fn=ebnm.fn,init.fn = init.fn)
    t3 = Sys.time()
    run_time_flash_greedy[iter] = difftime(t3,t2,units='secs')

    fit_flash = suppressWarnings(flash.backfit(fit_flash, verbose = verbose_flash,maxiter = maxiter_backfitting))
    t4 = Sys.time()
    run_time_flash_backfitting[iter] = difftime(t4,t3,units='secs')

    fit_flash = flash.nullcheck(fit_flash, verbose = verbose_flash)
    t5 = Sys.time()
    run_time_flash_nullcheck[iter] = difftime(t5,t4,units='secs')

    #print(fitted(fit_flash)[1:3,1:3])
    # fit_flash = flash.init(M, S = sqrt(sigma2), var.type = NULL, S.dim = S.dim)%>%
    #   flash.init.factors(init = fit_flash) %>%
    #   flash.add.greedy(Kmax = Kmax,verbose = verbose_flash) %>%
    #   flash.backfit(verbose = verbose_flash,maxiter = maxiter_backfitting) %>%
    #   flash.nullcheck(verbose = verbose_flash)

    K_trace[iter+1] = fit_flash$n.factors

    # fit_flash = flash.init(M, S = sqrt(sigma2), var.type = NULL, S.dim = S.dim)%>%
    #   flash.add.greedy(Kmax = Kmax,verbose = verbose_flash)


    KL_LF = sum(ff.KL(fit_flash$flash.fit,1)) + sum(ff.KL(fit_flash$flash.fit,2))
    #print(paste('After flash,elbo is',calc_split_PMF_obj_flashier(Y,S,sigma2,M,V,fit_flash,KL_LF,const,var_type)))




    # check convergence
    obj[iter + 1] = calc_split_PMF_obj_flashier(Y,S,sigma2,M,V,fit_flash,KL_LF,const,var_type)
    if((obj[iter+1] - obj[iter])/num_points < conv_tol){
      if((obj[iter+1] - obj[iter])<0){
        warning('An iteration decreases ELBO')
      }
      break
    }

    if(verbose){
      if(iter%%printevery==0){
        print(paste('At iter ',iter, ', ELBO=',round(obj[iter+1],log10(1/conv_tol)),sep = ''))
      }
    }

    # save fitted values because the total running time could be very long

    if(iter%%save_fit_every==0){
      saveRDS(list(fit_flash=fit_flash,
                   elbo=obj[length(obj)],
                   K_trace=K_trace,
                   elbo_trace=obj,
                   sigma2 = sigma2,
                   run_time = difftime(Sys.time(),start_time,units='auto'),
                   #M=M,V=V,
                   #init_val=init_val,
                   run_time_break_down = list(run_time_vga_init = run_time_vga_init,
                                              run_time_vga = run_time_vga,
                                              run_time_flash_init = run_time_flash_init,
                                              run_time_flash_init_factor = run_time_flash_init_factor,
                                              run_time_flash_greedy = run_time_flash_greedy,
                                              run_time_flash_backfitting = run_time_flash_backfitting,
                                              run_time_flash_nullcheck = run_time_flash_nullcheck)),
              file=paste(save_fit_path,save_fit_name,'_fit_iter',iter,'.rds',sep=''))
    }


  }

  end_time = Sys.time()
  return(list(fit_flash=fit_flash,
              elbo=obj[length(obj)],
              K_trace=K_trace,
              elbo_trace=obj,
              sigma2 = sigma2,
              run_time = difftime(end_time,start_time,units='auto'),
              #M=M,V=V,
              init_val=init_val,
              run_time_break_down = list(run_time_vga_init = run_time_vga_init,
                                         run_time_vga = run_time_vga,
                                         run_time_flash_init = run_time_flash_init,
                                         run_time_flash_init_factor = run_time_flash_init_factor,
                                         run_time_flash_greedy = run_time_flash_greedy,
                                         run_time_flash_backfitting = run_time_flash_backfitting,
                                         run_time_flash_nullcheck = run_time_flash_nullcheck)))
}


#'@title calc splitting PMF objective function.
calc_split_PMF_obj_flashier = function(Y,S,sigma2,M,V,fit_flash,KL_LF,const,var_type){
  R2 = fit_flash$flash.fit$R2
  n = nrow(Y)
  p = ncol(Y)
  if(var_type=='by_row'){
    sv = rowSums(V)
    ss = p
  }else if(var_type=='by_col'){
    sv = colSums(V)
    ss = n
  }else if(var_type=='constant'){
    sv = sum(V)
    ss = n*p
  }else{
    stop('Non-supported var type')
  }
  val = sum(Y*M - S*exp(M+V/2)   + log(2*pi*V)/2 + 0.5 ) - sum(ss*log(2*pi*sigma2)/2)- sum(sv/2/sigma2) - sum(R2/2/sigma2) + const+ KL_LF
  return(val)
}


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


#'@title a matrix version of the solver
#'@importFrom vebpm vga_pois_solver
vga_pois_solver_mat = function(init_Val,X,S,Beta,Sigma2,maxiter=1000,tol=1e-8){

  n = nrow(X)
  p = ncol(X)
  # transform them to vector
  x = as.vector(X)
  s = as.vector(S)
  beta = as.vector(Beta)
  sigma2 = as.vector(Sigma2)
  init_val = as.vector(init_Val)

  res = vga_pois_solver(init_val,x,s,beta,sigma2,maxiter=maxiter,tol=tol,method='newton')

  return(list(M = matrix(res$m,nrow=n,ncol=p,byrow = F),V = matrix(res$v,nrow=n,ncol=p,byrow = F)))

}
