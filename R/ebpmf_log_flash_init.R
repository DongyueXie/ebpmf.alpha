#'@title the init fit of flash in ebpmf log function
#'@importFrom ebnm ebnm_normal
ebpmf_log_flash_init = function(M,sigma2,l0,f0,ones_n,ones_p,loadings_sign,factors_sign,ebnm.fn,ebnm.fn.offset,
                                S.dim,verbose_flash,fix_l0,fix_f0,Kmax,add_greedy_extrapolate,maxiter_backfitting,
                                backfit_extrapolate,backfit_warmstart,
                                init.fn.flash,no_backfit_kset,n_refit_max){

  n = nrow(M)
  p = ncol(M)

  ebnm.fixed.l0 = function(x,s,g_init,fix_g,output){
    return(list(posterior=list(mean=l0,second_moment = l0^2),
                fitted_g = NULL,
                log_likelihood=sum(dnorm(x,l0,s,log=T))))
  }

  ebnm.fixed.f0 = function(x,s,g_init,fix_g,output){
    return(list(posterior=list(mean=f0,second_moment = f0^2),
                fitted_g = NULL,
                log_likelihood=sum(dnorm(x,f0,s,log=T))))
  }

  fit_flash = flash.init(M,S=sqrt(sigma2),var.type = NULL, S.dim=S.dim)%>%
              flash.set.verbose(verbose_flash)

  if(fix_l0){
    fit_flash = flash.init.factors(fit_flash,list(l0, ones_p),ebnm.fn = ebnm.fixed.l0) %>%
                flash.fix.factors(kset = 1, mode = 2)
  }else{
    fit_flash = flash.init.factors(fit_flash,list(l0, ones_p),ebnm.fn = ebnm.fn.offset) %>%
      flash.fix.factors(kset = 1, mode = 2)
  }

  if(fix_f0){
    fit_flash = flash.init.factors(fit_flash,list(ones_n, f0),ebnm.fn = ebnm.fixed.f0) %>%
      flash.fix.factors(kset = 2, mode = 1)
  }else{
    fit_flash = flash.init.factors(fit_flash,list(ones_n, f0),ebnm.fn = ebnm.fn.offset) %>%
      flash.fix.factors(kset = 2, mode = 1)
  }

  fit_flash = flash.add.greedy(fit_flash, Kmax = Kmax,ebnm.fn = ebnm.fn,init.fn=init.fn.flash,extrapolate = add_greedy_extrapolate)
  kset_backfit = (1:fit_flash$n.factors)[!(1:fit_flash$n.factors)%in%no_backfit_kset]
  fit_flash = suppressWarnings(flash.backfit(fit_flash,kset = kset_backfit,maxiter = maxiter_backfitting,extrapolate=backfit_extrapolate,warmstart = backfit_warmstart)%>%
                                 flash.nullcheck(kset=kset_backfit))

  n_refit = 0
  while(fit_flash$n.factors<=2&n_refit<=n_refit_max){
    n_refit = n_refit + 1
    cat(paste('No new structure found yet. Re-trying...',n_refit))
    cat('\n')
    init.fn.flash = function(f){init.fn.default(f, dim.signs = c(loadings_sign, factors_sign),seed = n_refit)}
    fit_flash = flash.add.greedy(fit_flash, Kmax = Kmax,ebnm.fn = ebnm.fn,init.fn=init.fn.flash,extrapolate = add_greedy_extrapolate)
    if(n_refit==n_refit_max){
      fit_flash = suppressWarnings(flash.backfit(fit_flash,kset = kset_backfit,maxiter = maxiter_backfitting,extrapolate=backfit_extrapolate,warmstart = backfit_warmstart))
    }else{
      fit_flash = suppressWarnings(flash.backfit(fit_flash,kset = kset_backfit,maxiter = maxiter_backfitting,extrapolate=backfit_extrapolate,warmstart = backfit_warmstart)%>%
                                     flash.nullcheck(kset=kset_backfit))
    }
  }

  if(fit_flash$n.factors<=2){
    warning('No new structure found in initialization.')
  }

  return(fit_flash)
}


