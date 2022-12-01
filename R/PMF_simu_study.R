#'@title a simulation study function
#'@importFrom parallel mclapply
#'@export
simu_study_PMF = function(simdata,n_cores = 1,method_list=c('flash','splitting'),Kmax=10,var_type='by_row'){
  n_simu = simdata$params$n_simu
  n_method = length(method_list)
  res = mclapply(1:n_simu,function(i){


    fitted_model <- vector("list", n_method)
    names(fitted_model) <- method_list
    S0 = tcrossprod(simdata$L0[,i],simdata$F0[,i])
    fitted_model$flash = try(flash(log(1+simdata$Y[,,i]/S0),Kmax = Kmax,var_type=var_type,verbose = F))
    fitted_model$splitting = try(splitting_PMF(simdata$Y[,,i],S0,Kmax=Kmax,var_type=var_type))

    mse_log = NULL
    k_hat = NULL
    if(class(fitted_model$flash)!='try-error'){
      mse_log = c(mse_log,rmse(fitted_model$flash$fitted_values,tcrossprod(simdata$Loading[,,i],simdata$Factor[,,i])))
      k_hat = c(k_hat,fitted_model$flash$nfactors)
    }else{
      mse_log = c(mse_log,NA)
      k_hat = c(k_hat,NA)
    }
    if(class(fitted_model$splitting)!='try-error'){
      mse_log = c(mse_log,rmse(fitted_model$splitting$fit_flash$fitted_values,tcrossprod(simdata$Loading[,,i],simdata$Factor[,,i])))
      k_hat = c(k_hat,fitted_model$splitting$fit_flash$nfactors)
    }else{
      mse_log = c(mse_log,NA)
      k_hat = c(k_hat,NA)
    }
    names(mse_log) = c('flash','splitting')
    names(k_hat) = c('flash','splitting')

    return(list(fitted_model = fitted_model,
                rmse_log=mse_log,
                k_hat=k_hat))

  },mc.cores = n_cores)
  return(list(sim_data=simdata,output = res))
}
