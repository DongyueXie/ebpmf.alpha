

#compare methods like nmf, stm-bmsm, stm-smashgen, stm-smashgen robust, HALS-wavelet, HALS-runmed

gene_splicing_study = function(X,name,K,nreps = 3,seed=12345){

  set.seed(seed)

  nmf_loss = Inf
  stm_loss = Inf
  stm_nugget_loss = Inf
  stm_nugget_robust_loss = Inf
  hals_wavelet_loss = Inf
  hals_runmed_loss = Inf

  for (reps in 1:nreps) {

    print(reps)

    fit_NMF = nnmf(X,k=K,method = 'scd',loss='mkl',verbose = F)
    if(min(fit_NMF$mkl)<nmf_loss){
      nmf_loss = min(fit_NMF$mkl)
      save(fit_NMF,file = paste('~/SMF/data/',name,'_NMF_mkl_scd_K',K,'.RData',sep = ''))
    }

    fit_stm = stm(X,K,init = list(L_init=fit_NMF$W,F_init = fit_NMF$H),
                  return_all = FALSE,tol=1e-2,maxiter=50,printevery = 1e5)
    if(fit_stm$KL[length(fit_stm$KL)]<stm_loss){
      stm_loss=fit_stm$KL[length(fit_stm$KL)]
      save(fit_stm,file=paste('~/SMF/data/',name,'_stm_bmsm_K',K,'.RData',sep = ''))
    }

    fit_stm_nugget_robust = stm(X,K,init = list(L_init=fit_NMF$W,F_init = fit_NMF$H),
                         return_all = FALSE,tol=1e-2,nugget = TRUE,maxiter=50,printevery = 1e5)
    if(fit_stm_nugget_robust$KL[length(fit_stm_nugget_robust$KL)]<stm_nugget_robust_loss){
      stm_nugget_robust_loss=fit_stm_nugget_robust$KL[length(fit_stm_nugget_robust$KL)]
      save(fit_stm_nugget_robust,file=paste('~/SMF/data/',name,'_stm_nugget_robust_K',K,'.RData',sep = ''))
    }

    fit_stm_nugget = stm(X,K,init = list(L_init=fit_NMF$W,F_init = fit_NMF$H),
                                return_all = FALSE,tol=1e-2,nugget = TRUE,maxiter=50,printevery = 1e5,
                         nug_control_f = list(robust=F))
    if(fit_stm_nugget$KL[length(fit_stm_nugget$KL)]<stm_nugget_loss){
      stm_nugget_loss=fit_stm_nugget$KL[length(fit_stm_nugget$KL)]
      save(fit_stm_nugget,file=paste('~/SMF/data/',name,'_stm_nugget_K',K,'.RData',sep = ''))
    }

    fit_hals_wavelet = NMF_HALS(X,K,smooth_method = 'wavelet',printevery = 1e5)
    if(fit_hals_wavelet$loss[length(fit_hals_wavelet$loss)]<hals_wavelet_loss){
      hals_wavelet_loss = fit_hals_wavelet$loss[length(fit_hals_wavelet$loss)]
      save(fit_hals_wavelet,file=paste('~/SMF/data/',name,'_hals_wavelet_K',K,'.RData',sep = ''))
    }

    fit_hals_runmed = NMF_HALS(X,K,smooth_method = 'runmed',printevery = 1e5)

    if(fit_hals_runmed$loss[length(fit_hals_runmed$loss)]<hals_runmed_loss){
      hals_runmed_loss = fit_hals_runmed$loss[length(fit_hals_runmed$loss)]
      save(fit_hals_runmed,file=paste('~/SMF/data/',name,'_hals_runmed_K',K,'.RData',sep = ''))
    }
  }
}

