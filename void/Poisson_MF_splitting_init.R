#'@title this function initializes the splitting PMF
#'@description sigma2, M needs initialization.
#'

splitting_PMF_init = function(Y,S,sigma2,
                      var_type,
                      M_init,
                      verbose,
                      n_cores,
                      init_tol,
                      printevery){
  n = nrow(Y)
  p = ncol(Y)
  # we need a. M; b. sigma2
  # case 1: both are given in input.
  if(!is.null(sigma2) & !is.null(M_init)){
    return(list(sigma2_init=sigma2,M_init=M_init,V_init = NULL))
  }

  # then we fit poisGG with mean = 0, var =sigma2
  if(is.null(M_init)){
    if(verbose){
      cat('Initializing M...')
    }
    # pre-estimate sigma2, assuming LF = 0?.
    if(var_type=='constant'){
      if(verbose){
        cat('Solving VGA...')
      }
      init_val = suppressWarnings(pois_mean_GG(as.vector(Y),as.vector(S),prior_mean = 0,prior_var = sigma2,tol=init_tol))
      M0 = matrix(init_val$posterior$mean_log,nrow=n,ncol=p)
      V = matrix(init_val$posterior$var_log,nrow=n,ncol=p)
      sigma2_init = init_val$fitted_g$var
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
        fit = suppressWarnings(pois_mean_GG(Y[i,],S[i,],prior_mean = 0,prior_var = sigma2[i],tol=init_tol))
        return(list(sigma2 = fit$fitted_g$var,mean_log = fit$posterior$mean_log,var_log = fit$posterior$var_log))
      },mc.cores = n_cores)
      sigma2_init = unlist(lapply(init_val,function(fit){fit$sigma2}))
      M0 = do.call(rbind,lapply(init_val,function(fit){fit$mean_log}))
      V = do.call(rbind,lapply(init_val,function(fit){fit$var_log}))
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
        fit = suppressWarnings(pois_mean_GG(Y[,i],S[,i],prior_mean = 0,prior_var = sigma2[i],tol=init_tol))
        return(list(sigma2 = fit$fitted_g$var,mean_log = fit$posterior$mean_log,var_log = fit$posterior$var_log))
      },mc.cores = n_cores)
      sigma2_init = unlist(lapply(init_val,function(fit){fit$sigma2}))
      M0 = do.call(cbind,lapply(init_val,function(fit){fit$mean_log}))
      V = do.call(cbind,lapply(init_val,function(fit){fit$var_log}))
    }
    return(list(sigma2_init=sigma2_init,M_init=M0,V_init = V))
  }
  # case 3: M_init is given, sigma2 is unknown, this case is not interesting and very unrealistic even for development purpose.

}

