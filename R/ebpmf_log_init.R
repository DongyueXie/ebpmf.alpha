#'@title this function initializes the ebpmf with log link
#'@description sigma2, M needs initialization.
#'

ebpmf_log_init = function(Y,l0,f0,sigma2,
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
    return(list(sigma2_init=sigma2,M_init=M_init))
  }

  # then we fit poisGG with mean = 0, var =sigma2
  if(is.null(M_init)){
    if(verbose){
      cat('Initializing M...')
    }
    # pre-estimate sigma2, assuming LF = 0?.
    if(var_type=='constant'){
      if(verbose){
        cat('Solving VGA constant...For large matrix this may require large memory usage')
      }
      init_val = suppressWarnings(pois_mean_GG(as.vector(Y),as.vector(tcrossprod(l0,f0)),prior_mean = 0,prior_var = sigma2,tol=init_tol))
      M = matrix(init_val$posterior$mean_log,nrow=n,ncol=p)
      sigma2_init = init_val$fitted_g$var
      rm(init_val)
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
        fit = suppressWarnings(pois_mean_GG(Y[i,],l0[i]*f0,prior_mean = 0,prior_var = sigma2[i],tol=init_tol))
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
          fit = suppressWarnings(pois_mean_GG(Y[,i],l0*f0[i],prior_mean = 0,prior_var = sigma2[i],tol=init_tol))
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
          fit = suppressWarnings(pois_mean_GG(Y[,j],l0*f0[j],prior_mean = 0,prior_var = sigma2[j],tol=init_tol))
          M[,j] = fit$posterior$mean_log
          sigma2_init[j] = fit$fitted_g$var
        }
        rm(fit)
      }
    }
    return(list(sigma2_init=sigma2_init,M_init=M))
  }
  ## case 3: M_init is given, sigma2 is unknown, this case is useful when using glmpca for initialization.
  ## In this case we will run flash to init sigma2
  return(list(sigma2_init=NULL,M_init=M_init))
  # if(!is.null(M_init)){
  #   if(verbose){
  #     cat('Initializing sigma2...')
  #   }
  #
  #   if(var_type=='constant'){
  #     if(verbose){
  #       cat('Solving VGA constant...For large matrix this may require large memory usage')
  #     }
  #     init_val = suppressWarnings(pois_mean_GG(as.vector(Y),as.vector(tcrossprod(l0,f0)),prior_mean = as.vector(M_init),prior_var = sigma2,tol=init_tol))
  #     sigma2_init = init_val$fitted_g$var
  #     rm(init_val)
  #   }
  #   if(var_type=='by_row'){
  #     if(verbose){
  #       cat('Solving VGA for row 1...')
  #     }
  #     init_val = mclapply(1:n,function(i){
  #       if(verbose){
  #         if(i%%printevery==0){
  #           cat(paste(i,'...'))
  #         }
  #       }
  #       fit = suppressWarnings(pois_mean_GG(Y[i,],l0[i]*f0,prior_mean = M_init[i,],prior_var = sigma2[i],tol=init_tol))
  #       return(fit$fitted_g$var)
  #     },mc.cores = n_cores)
  #     sigma2_init = unlist(init_val)
  #     rm(init_val)
  #   }
  #   if(var_type=='by_col'){
  #     if(verbose){
  #       cat('Solving VGA for column 1...')
  #     }
  #
  #     if(n_cores>1){
  #       init_val = mclapply(1:p,function(i){
  #         if(verbose){
  #           if(i%%printevery==0){
  #             cat(paste(i,'...'))
  #           }
  #         }
  #         fit = suppressWarnings(pois_mean_GG(Y[,i],l0*f0[i],prior_mean = M_init[,i],prior_var = sigma2[i],tol=init_tol))
  #         return(fit$fitted_g$var)
  #       },mc.cores = n_cores)
  #       sigma2_init = unlist(init_val)
  #       rm(init_val)
  #     }else{
  #       M = matrix(nrow=n,ncol=p)
  #       sigma2_init = rep(0,p)
  #       for(j in 1:p){
  #         if(verbose){
  #           if(j%%printevery==0){
  #             cat(paste(j,'...'))
  #           }
  #         }
  #         fit = suppressWarnings(pois_mean_GG(Y[,j],l0*f0[j],prior_mean = M_init[,j],prior_var = sigma2[j],tol=init_tol))
  #         sigma2_init[j] = fit$fitted_g$var
  #       }
  #       rm(fit)
  #     }
  #   }
  #   return(list(sigma2_init=sigma2_init,M_init=M_init))
  # }
}

