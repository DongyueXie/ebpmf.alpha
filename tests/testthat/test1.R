set.seed(12345)
N = 100
p = 50
K = 2
Ftrue = matrix(0,nrow=p,ncol=K)
Ftrue[1:20,1] = 1
Ftrue[21:40,2] = 1
Ltrue = matrix(rnorm(N*K), ncol=K)
# test
Lambda = exp(tcrossprod(Ltrue,Ftrue))
Y = matrix(rpois(N*p,Lambda),nrow=N,ncol=p)
fit = splitting_PMF_flashier(Y,verbose=TRUE)
plot(fitted(fit$fit_flash),tcrossprod(Ltrue,Ftrue),col='grey80')
abline(a=0,b=1)
# test nonegative loading option
Lambda = exp(tcrossprod(abs(Ltrue),Ftrue))
Y = matrix(rpois(N*p,Lambda),nrow=N,ncol=p)
fit = splitting_PMF_flashier(Y,verbose=TRUE,
                             ebnm.fn = c(ebnm::ebnm_point_exponential, ebnm::ebnm_point_normal),
                             loadings_sign = 1,maxiter = 100)
plot(fitted(fit$fit_flash),tcrossprod(abs(Ltrue),Ftrue),col='grey80')
abline(a=0,b=1)
plot(fit$fit_flash$F.pm[,1],type='l')
plot(fit$fit_flash$F.pm[,2],type='l')

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
