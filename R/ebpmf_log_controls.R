ebpmf_log_init_control_default = function(){
  return(list(n_cores=1,
              init_tol = 1e-5,
              verbose = TRUE,
              printevery = 100,
              single_gene_ebpm = TRUE,
              sigma2_init = NULL,
              M_init = NULL
              ))
}


ebpmf_log_general_control_default = function(){
  return(list(batch_size = Inf,
              maxiter=50,
              conv_tol=1e-5,
              printevery=10,
              verbose=TRUE,
              save_init_val = FALSE,
              save_latent_M = FALSE,
              garbage_collection_every = 10,
              save_fit_every = Inf,
              save_fit_path = NULL,
              save_fit_name = NULL))
}

ebpmf_log_vga_control_default = function(){
  return(list(maxiter_vga = 10,
         vga_tol = 1e-5))
}

ebpmf_log_sigma2_control_default = function(){
  return(list(est_sigma2 = TRUE,
              a0 = 1,
              b0 = 1,
              cap_var_mean_ratio = 0,
              return_sigma2_trace = FALSE))
}

#'This function returns init.fn.flash and no_backfit_kset
flash_extra_control = function(loadings_sign,factors_sign,fix_l0,fix_f0){
  if(loadings_sign==0&factors_sign==0){
    # this is faster than the default init method in flash
    init.fn.flash = function(f){init.fn.irlba(f)}
  }else{
    init.fn.flash = function(f){init.fn.default(f, dim.signs = c(loadings_sign, factors_sign))}
  }
  no_backfit_kset = NULL
  if(fix_l0){
    no_backfit_kset = c(no_backfit_kset,1)
  }
  if(fix_f0){
    no_backfit_kset = c(no_backfit_kset,2)
  }
  return(list(init.fn.flash=init.fn.flash,no_backfit_kset=no_backfit_kset))
}

#'@title default flash parameters
#'@importFrom ebnm ebnm_point_normal
ebpmf_log_flash_control_default = function(){
  return(list(ebnm.fn = ebnm_point_normal,
              ebnm.fn.offset = ebnm_normal,
              loadings_sign = 0,
              factors_sign = 0,
              fix_l0=TRUE,
              fix_f0=FALSE,
              Kmax=30,
              add_greedy_Kmax = 1,
              add_greedy_warmstart = TRUE,
              add_greedy_extrapolate = FALSE,
              add_greedy_every = 1,
              backfit_extrapolate = TRUE,
              backfit_warmstart = TRUE,
              maxiter_backfitting = 1,
              verbose_flash=0
  ))
}
