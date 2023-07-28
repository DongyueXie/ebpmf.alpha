sla_2000 <- readRDS("~/Downloads/sla_2000.rds")
Y = sla_2000$data

l0 = log(cbind(rowMeans(Y)))
f0 = log(cbind(colSums(Y)/sum(exp(l0))))
devtools::load_all()
n_cores = 5
ebpm_init = FALSE
init_tol = 1e-5
init_maxiter = 100

p = ncol(Y)
n = nrow(Y)

init_val = mclapply(1:p,function(j){
  if(T){
    if(j%%10==0){
      cat(paste(j,'...'))
    }
  }
  if(ebpm_init){
    fit = ash_pois(Y[,j],scale = drop(exp(l0+f0[j])),link='identity',mode=1,method='shrink')
    #init_var_vga = mean(log(fit$result$PosteriorMean)^2)
    init_var_vga = NULL
    fit = suppressWarnings(ebpm_normal(Y[,j],
                                       g_init = list(mean=drop(l0+f0[j]),var=init_var_vga),
                                       q_init = list(m_init=log(fit$result$PosteriorMean) + drop(l0+f0[j]),v_init = init_var_vga),
                                       fix_g = c(TRUE,FALSE),tol=init_tol,conv_type = conv_type,
                                       maxiter = init_maxiter))
  }else{
    fit = suppressWarnings(ebpm_normal(Y[,j],
                                       g_init = list(mean=drop(l0+f0[j]),var=NULL),
                                       fix_g = c(TRUE,FALSE),tol=init_tol,
                                       maxiter = init_maxiter))
  }

  return(list(sigma2 = fit$fitted_g$var,mean_log = fit$posterior$mean_log))
},mc.cores = n_cores)
sigma2_init = unlist(lapply(init_val,function(fit){fit$sigma2}))
M = do.call(cbind,lapply(init_val,function(fit){fit$mean_log}))

hist(sigma2_init,breaks = 100)
alpha = l0%*%t(rep(1,p)) + rep(1,n)%*%t(f0)
plot(colSums(Y),sigma2_init)
plot(colSums(exp(alpha)),sigma2_init)
plot((colSums(Y/exp(alpha))),sigma2_init)

plot(M[,2],l0+f0[2])
plot(M[,2],Y[,2])
abline(a=0,b=1)

which(Y[,2]==3)
log(Y[328,2]/exp(l0[328]+f0[2]))
ebpm_normal(Y[328,2],s=exp(drop(l0[328]+f0[2])),
            g_init = list(mean=0,var=2.716657^2),
            fix_g = c(T,T),tol=init_tol,
            maxiter = init_maxiter)

temp = ebpm_normal(Y[,2],s=exp(drop(l0+f0[2])),
                   g_init = list(mean=0,var=100),
                   fix_g = c(T,T),tol=init_tol,
                   maxiter = init_maxiter)
plot(temp$posterior$mean_log)

plot(M[,1],Y[,1])



# fit EBMF on M
fit = flash(M-alpha,ebnm.fn = ebnm::ebnm_point_exponential,var.type = 2,greedy.Kmax = 6,backfit = T)
colnames(sla_2000$data)[order(fit$F.pm[,1],decreasing = T)[1:10]]
colnames(sla_2000$data)[order(fit$F.pm[,2],decreasing = T)[1:10]]
colnames(sla_2000$data)[order(fit$F.pm[,3],decreasing = T)[1:10]]
colnames(sla_2000$data)[order(fit$F.pm[,4],decreasing = T)[1:10]]
colnames(sla_2000$data)[order(fit$F.pm[,5],decreasing = T)[1:10]]
colnames(sla_2000$data)[order(fit$F.pm[,6],decreasing = T)[1:10]]

fit = flash(M-alpha,ebnm.fn = ebnm::ebnm_point_exponential,var.type = 0,S=sqrt(sigma2_init),greedy.Kmax = 6,backfit = T)



fit0 = ebpmf_log(sla_2000$data,flash_control=list(Kmax=6,
                                              ebnm.fn=c(ebnm::ebnm_point_exponential,ebnm::ebnm_point_exponential),
                                              loadings_sign=1,
                                              factors_sign = 1),
                var_type = 'by_col',
                init_control = list(init_tol=1e-5,deal_with_no_init_factor='flash_dryrun',n_cores=5,M_init=fit0$init_val$M_init,sigma2_init = fit0$init_val$sigma2_init),
                sigma2_control = list(return_sigma2_trace=T),
                general_control = list(maxiter=10,conv_tol=1e-5,save_init_val=TRUE,save_latent_M=T))

fit = flash(fit0$init_val$M_init-alpha,ebnm.fn = ebnm::ebnm_point_exponential,var.type = 2,greedy.Kmax = 6,backfit = T)

#################################################

init_val = mclapply(1:p,function(j){
  if(T){
    if(j%%10==0){
      cat(paste(j,'...'))
    }
  }
  if(ebpm_init){
    fit = ash_pois(Y[,j],scale = drop(exp(l0+f0[j])),link='identity',mode=1,method='shrink')
    #init_var_vga = mean(log(fit$result$PosteriorMean)^2)
    init_var_vga = NULL
    fit = suppressWarnings(ebpm_normal(Y[,j],
                                       g_init = list(mean=drop(l0+f0[j]),var=init_var_vga),
                                       q_init = list(m_init=log(fit$result$PosteriorMean) + drop(l0+f0[j]),v_init = init_var_vga),
                                       fix_g = c(TRUE,FALSE),tol=init_tol,conv_type = conv_type,
                                       maxiter = init_maxiter))
  }else{
    fit = suppressWarnings(ebpm_normal(Y[,j],s = exp(drop(l0+f0[j])),
                                       g_init = list(mean=0,var=max(log(Y[,j]/exp(l0+f0[j])))^2),
                                       fix_g = c(F,TRUE),tol=init_tol,
                                       maxiter = init_maxiter))
  }

  return(list(sigma2 = fit$fitted_g$var,mean_log = fit$posterior$mean_log))
},mc.cores = n_cores)
sigma2_init = unlist(lapply(init_val,function(fit){fit$sigma2}))
M = do.call(cbind,lapply(init_val,function(fit){fit$mean_log}))

fit = flash(M,ebnm.fn = ebnm::ebnm_point_exponential,var.type = 2,greedy.Kmax = 6,backfit = T)

fit2 = NNLM::nnmf(M,k=6,loss='mse')
plot(fit2$H[6,])

for(k in 1:6){
  print(colnames(sla_2000$data)[order(fit2$H[k,],decreasing = T)[1:10]])
}

Y_tilde = log_for_ebmf(sla_2000$data)
Y_tilde2  = Y_tilde - rep(1,nrow(Y_tilde))%*%t(colMeans(Y_tilde))

plot(Y_tilde2[,3])
plot(M[,3])
plot(Y[,3])

Y_tilde0 = log(1e-5+Y/exp(alpha))
fit = flash(Y_tilde0,ebnm.fn = ebnm::ebnm_point_exponential,var.type = 2,greedy.Kmax = 6,backfit = T)

for(i in 1:fit$n.factors){
  print(colnames(sla_2000$data)[order(fit$F.pm[,i],decreasing = T)[1:10]])
}


D = Y*log((Y+0.0001)/exp(alpha)) - Y + exp(alpha)
fit = flash(D,ebnm.fn = ebnm::ebnm_point_exponential,var.type = 2,greedy.Kmax = 6,backfit = T)

for(i in 1:fit$n.factors){
  print(colnames(sla_2000$data)[order(fit$F.pm[,i],decreasing = T)[1:10]])
}

###########
fit_ebpmf = ebpmf_log(sla_2000$data,
                      flash_control=list(backfit_extrapolate=T,backfit_warmstart=T,
                                         ebnm.fn = c(ebnm::ebnm_point_exponential, ebnm::ebnm_point_exponential),
                                                       loadings_sign = 1,factors_sign=1,Kmax=10),
                      init_control = list(n_cores=5,flash_est_sigma2=F,log_init_for_non0y=T))


for(i in 3:fit_ebpmf$fit_flash$n_factors){
  cat(colnames(sla_2000$data)[order(fit_ebpmf$fit_flash$F_pm[,i],decreasing = T)[1:10]]);cat('\n')
}
plot(fit_ebpmf$elbo_trace)
plot(fit_ebpmf$fit_flash$F_pm[,4])
hist(fit_ebpmf$sigma2,breaks = 100)

fit_ebpmf2 = ebpmf_log(sla_2000$data,
                      flash_control=list(backfit_extrapolate=T,backfit_warmstart=T,
                                         ebnm.fn = c(ebnm::ebnm_point_exponential, ebnm::ebnm_point_exponential),
                                         loadings_sign = 1,factors_sign=1,Kmax=10),
                      init_control = list(n_cores=5,flash_est_sigma2=T,log_init_for_non0y=T))

fit_ebpmf3 = ebpmf_log(sla_2000$data,
                       flash_control=list(backfit_extrapolate=T,backfit_warmstart=T,
                                          ebnm.fn = c(ebnm::ebnm_point_exponential, ebnm::ebnm_point_exponential),
                                          loadings_sign = 1,factors_sign=1,Kmax=10),
                       init_control = list(n_cores=5,flash_est_sigma2=T,log_init_for_non0y=F))

# fit_ebpmf4 = ebpmf_log(sla_2000$data,
#                        flash_control=list(backfit_extrapolate=F,backfit_warmstart=F,
#                                           ebnm.fn = c(ebnm::ebnm_point_exponential, ebnm::ebnm_point_exponential),
#                                           loadings_sign = 1,factors_sign=1,Kmax=10),
#                        init_control = list(n_cores=5,flash_est_sigma2=F,log_init_for_non0y=F))
#
#

# fit_ebpmf2 = ebpmf_log(sla_2000$data,
#                       flash_control=list(backfit_extrapolate=FALSE,backfit_warmstart=F,ebnm.fn = c(ebnm::ebnm_point_exponential, ebnm::ebnm_point_normal),
#                                          loadings_sign = 1,factors_sign=0),init_control = list(n_cores=5))
#
