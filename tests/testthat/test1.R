set.seed(12345)
N = 100
p = 50
K = 2
Ftrue = matrix(0,nrow=p,ncol=K)
Ftrue[1:20,1] = 1
Ftrue[21:40,2] = 1
Ltrue = matrix(rnorm(N*K), ncol=K)

l0 = runif(N,1,2)
f0 = runif(p,1,2)
S0 = tcrossprod(l0,f0)

Lambda = S0*exp(tcrossprod(Ltrue,Ftrue))

Y = matrix(rpois(N*p,Lambda),nrow=N,ncol=p)

fit = splitting_PMF(Y,S0)



simdata= sim_data_log(n=100,p=100,K=3,n_simu = 2)
res = simu_study_PMF(simdata)
