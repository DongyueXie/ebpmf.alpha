load('/project2/mstephens/pcarbo/git/single-cell-topics/data/pbmc_purified.RData')
library(Matrix)
l0 = c(rowSums(counts))/sqrt(sum(counts))
f0 = c(colSums(counts))/sqrt(sum(counts))

library(stm)
library(peakRAM)
peakRAM(fit <- splitting_PMF_flashier_low_memory(counts,l0,f0,printevery = 1,
                                                 verbose = TRUE,maxiter = 500,
                                                 vga_tol=1e-3,
                                                 init_tol=1e-3,
                                                 conv_tol = 1e-3,
                                                 n_cores = 1))
saveRDS(fit,'pbmc_purified.rds')

M = as.matrix(log(0.5+counts))
gc()
fit_flash = flash(M,var.type=2,greedy.Kmax =3)
EL = fit_flash$L.pm
EF = fit_flash$F.pm
sigma2 = fit_flash$residuals.sd^2
gc()
idx = 1:10000
peakRAM::peakRAM(temp <- vga_pois_solver_mat_newton_low_memory(M[idx,],counts[idx,],l0[idx],f0,
                                                        EL[idx,],
                                                        EF,
                                                        sigma2,
                                                        var_type='by_col',
                                                        maxiter=1000,
                                                        tol=1e-5))

peakRAM::peakRAM(temp <- vga_pois_solver_mat_newton(M[idx,],counts[idx,],tcrossprod(l0[idx],f0),
                                                               tcrossprod(EL[idx,],EF),
                                                              matrix(sigma2,nrow=length(idx),ncol=ncol(M),byrow = TRUE),
                                                               maxiter=1000,
                                                               tol=1e-3,return_V = F))

m0 = M[idx,]
y0 = as.matrix(counts[idx,])
s0 = tcrossprod(l0[idx],f0)
beta0 = tcrossprod(EL[idx,],EF)
sigma20 = matrix(sigma2,nrow=length(idx),ncol=ncol(M),byrow = TRUE)
peakRAM::peakRAM(temp <- vga_pois_solver_mat_newton(m0,y0,s0,
                                                    beta0,
                                                    sigma20,
                                                    maxiter=5,
                                                    tol=1e-3,return_V = F))


##############################################
set.seed(12345)
N = 1000
p = 1000
K = 2
sigma2 = 0.1
Ftrue = matrix(0,nrow=p,ncol=K)
Ftrue[1:20,1] = 1
Ftrue[21:40,2] = 1
#Ftrue[41:60,3] = 1
Ltrue = matrix(rnorm(N*K), ncol=K)
# test
Lambda = exp(tcrossprod(Ltrue,Ftrue) + matrix(rnorm(N*p,0,sqrt(sigma2)),nrow=N))
Y = matrix(rpois(N*p,Lambda),nrow=N,ncol=p)
Y[Y<quantile(Y,0.98)] = 0
sum(Y==0)/prod(dim(Y))
library(peakRAM)
gc()
peakRAM(fit <- splitting_PMF_flashier(Y,verbose=TRUE,n_cores = 1,maxiter_vga = 3,maxiter = 3,printevery = 1,vga_tol = 0.1,init_tol = 0.1,Kmax_init=2))

gc()
Ys=Matrix(Y,sparse = T)
peakRAM(fit <- splitting_PMF_flashier_low_memory(Ys,verbose=TRUE,n_cores = 1,maxiter_vga = 3,printevery = 1,batch_size = 1000,maxiter = 3,vga_tol = 0.1,init_tol = 0.1,Kmax_init=2))
