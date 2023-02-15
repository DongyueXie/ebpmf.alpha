#'@title this function fits flash during ebpmf log iterations
#'@param fit_flash An object from running flash()
ebpmf_log_flash_update = function(fit_flash,sigma2,ones_n,ones_p,iter,loadings_sign,factors_sign,ebnm.fn,ebnm.fn.offset,
                                  S.dim,verbose_flash,fix_l0,fix_f0,Kmax,add_greedy_extrapolate,maxiter_backfitting,
                                  add_greedy_every,add_greedy_Kmax,add_greedy_warmstart,
                                  init.fn.flash,no_backfit_kset){

  ## create an init flash.fit obj for flash.init.factor
  ## use flash.fit for init can also init the posterior second moment
  ## Also init the prior g
  flash_fit_init = fit_flash$flash.fit[c('EF','EF2','g')]
  class(flash_fit_init) = c("flash.fit","list")


  l0 = flash_fit_init$EF[[1]][,1,drop=FALSE]
  f0 = flash_fit_init$EF[[2]][,2,drop=FALSE]
  flash_fit_init$EF[[1]] = flash_fit_init$EF[[1]][,-c(1,2),drop=FALSE]
  flash_fit_init$EF2[[1]] = flash_fit_init$EF2[[1]][,-c(1,2),drop=FALSE]
  flash_fit_init$EF[[2]] = flash_fit_init$EF[[2]][,-c(1,2),drop=FALSE]
  flash_fit_init$EF2[[2]] = flash_fit_init$EF2[[2]][,-c(1,2),drop=FALSE]

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

  fit_flash = flash.init(fit_flash$flash.fit$Y, S = sqrt(sigma2), var.type = NULL, S.dim = S.dim) %>%
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
  fit_flash = flash.init.factors(fit_flash,flash_fit_init,ebnm.fn = ebnm.fn)
  fit_flash$flash.fit$g = flash_fit_init$g
  if(iter%%add_greedy_every==0 & fit_flash$n.factors < Kmax){
    fit_flash = flash.add.greedy(fit_flash, Kmax = add_greedy_Kmax,
                                 ebnm.fn=ebnm.fn,init.fn = init.fn.flash,
                                 warmstart = add_greedy_warmstart,
                                 extrapolate = add_greedy_extrapolate)
  }

  fit_flash = suppressWarnings(flash.backfit(fit_flash, kset = (1:fit_flash$n.factors)[!(1:fit_flash$n.factors)%in%no_backfit_kset],maxiter = maxiter_backfitting)%>%
                                 flash.nullcheck())


  fit_flash


}
