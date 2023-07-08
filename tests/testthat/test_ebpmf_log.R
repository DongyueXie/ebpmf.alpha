library(peakRAM)
set.seed(12345)
N = 1000
p = 100
K = 3
sigma2 = 0
Ftrue = matrix(0,nrow=p,ncol=K)
Ftrue[1:20,1] = 1
Ftrue[21:40,2] = 2
Ftrue[41:60,3] = 3
Ltrue = matrix(rnorm(N*K), ncol=K)
# test
Lambda = exp(tcrossprod(Ltrue,Ftrue) + matrix(rnorm(N*p,0,sqrt(sigma2)),nrow=N))
Y = matrix(rpois(N*p,Lambda),nrow=N,ncol=p)
sum(Y!=0)/prod(dim(Y))

peakRAM(fit <- ebpmf_log(Y,l0=0,f0=0,flash_control=list(fix_f0=T),init_control=list(conv_type='sigma2abs')))

# peakRAM(fit <- splitting_PMF_flashier_low_memory(Y,verbose=TRUE,n_cores = 10,maxiter_vga = 2,printevery = 1,batch_size = 100))
plot(fit$K_trace)
plot(fitted(fit$fit_flash),tcrossprod(Ltrue,Ftrue),col='grey80')
abline(a=0,b=1)
fit$fit_flash$pve
for(k in 1:fit$fit_flash$n.factors){
  plot(fit$fit_flash$F.pm[,k],type='l')
}
# test nonegative loading option
set.seed(12345)
Lambda = exp(tcrossprod(abs(Ltrue),Ftrue)+ matrix(rnorm(N*p,0,sqrt(sigma2)),nrow=N))
Y = matrix(rpois(N*p,Lambda),nrow=N,ncol=p)

fit = ebpmf_log(Y,l0=0,f0=0,flash_control=list(ebnm.fn = c(ebnm::ebnm_point_exponential, ebnm::ebnm_point_normal),
                             loadings_sign = 1,fix_f0=T))
plot(fit$K_trace)
fit$fit_flash$pve
plot(fitted(fit$fit_flash),tcrossprod(abs(Ltrue),Ftrue),col='grey80')
abline(a=0,b=1)
for(k in 1:fit$fit_flash$n.factors){
  plot(fit$fit_flash$F.pm[,k],type='l')
}


# test nonegative loading and factor option
set.seed(12345)
Lambda = exp(tcrossprod(abs(Ltrue),abs(Ftrue))+ matrix(rnorm(N*p,0,sqrt(sigma2)),nrow=N))
Y = matrix(rpois(N*p,Lambda),nrow=N,ncol=p)
fit = ebpmf_log(Y,l0=0,f0=0,flash_control=list(ebnm.fn = c(ebnm::ebnm_point_exponential, ebnm::ebnm_point_exponential),
                                               loadings_sign = 1,factors_sign=1,fix_f0=T))
for(k in 1:fit$fit_flash$n.factors){
  plot(fit$fit_flash$F.pm[,k],type='l')
}
#########################
# use point_exponential prior always gives an error. I suspect it's related to the optimization method.
# I changed optmethod to trust and it did not show errors. But instead returned an increasing number of factors.
# unimodal_nonnegative prior works fine but very slow
temp_func = function(x,
                     s = 1,
                     mode = 0,
                     scale = "estimate",
                     g_init = NULL,
                     fix_g = FALSE,
                     output = output_default(),
                     optmethod = 'trust',
                     control = NULL){ebnm_point_exponential(
                       x,
                       s,
                       mode ,
                       scale ,
                       g_init ,
                       fix_g ,
                       output ,
                       optmethod ,
                       control
                     )}

fit = splitting_PMF_flashier(Y,verbose=TRUE,
                             ebnm.fn = c(ebnm::ebnm_point_exponential, ebnm::ebnm_point_exponential),
                             loadings_sign = 1,factors_sign = 1,maxiter = 100,n_cores = 10)
fit = splitting_PMF_flashier(Y,verbose=TRUE,
                             add_greedy_Kmax = 1,
                             maxiter_vga = 100,
                             vga_tol = 1e-8,
                             add_greedy_init = 'new_init',
                             ebnm.fn = c(ebnm::ebnm_unimodal_nonnegative, ebnm::ebnm_unimodal_nonnegative),
                             loadings_sign = 1,factors_sign = 1,maxiter = 100,n_cores = 10)
plot(fit$K_trace)
plot(fitted(fit$fit_flash),tcrossprod(abs(Ltrue),abs(Ftrue)),col='grey80')
abline(a=0,b=1)
plot(fit$fit_flash$F.pm[,1],type='l')
plot(fit$fit_flash$F.pm[,2],type='l')
plot(fit$fit_flash$F.pm[,3],type='l')
plot(fit$fit_flash$F.pm[,4],type='l')
#
# fit = ebpmf(Y,verbose=TRUE,
#             ebnm.fn = c(ebnm::ebnm_point_exponential, ebnm::ebnm_point_exponential),
#             loadings_sign = 1,factors_sign = 1,maxiter = 20)

# datax= sim_data_log(n=100,p=100,K=3,n_simu = 1)
# res = simu_study_PMF(datax)
#
# datax= sim_data_log(n=500,p=500,K=3,n_simu = 1)
# Y = datax$Y[,,1]
# S = tcrossprod(datax$L0[,1],datax$F0[,1])
#
# datax = sim_data_log_simple(500,400)
# Y = datax$Y
# S = tcrossprod(datax$l0,datax$f0)





myPCA <- PLNPCA(Y ~ 1, ranks = 1:3)
myPCA_ICL <- getBestModel(myPCA, "ICL")
plot(myPCA_ICL$rotation[,1])
plot(myPCA_ICL$rotation[,2])
plot(myPCA_ICL$rotation[,3])


