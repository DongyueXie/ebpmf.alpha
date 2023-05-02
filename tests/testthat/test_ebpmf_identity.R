set.seed(123)
n = 120
p = 300
K= 3
L = matrix(0, nrow=n, ncol=K)
FF = matrix(0, nrow=K, ncol=p)
L[1:(n/3),1] = 1
L[((n/3)+1):(2*n/3),2] = 1
L[((2*n/3)+1):n,3] = 1
L = L + matrix(runif(n*K,0,0.5),nrow=n)
FF[1,1:(p/3)] = 1+10
FF[2,((p/3)+1):(2*p/3)] = 1+10
FF[3,((2*p/3)+1):p] = 1+10
FF = FF + matrix(rnorm(p*K,0,1),ncol=p)
FF = pmax(FF,0)
lambda = L %*% FF
X = matrix(rpois(n=length(lambda),lambda),nrow=n)

fit0 = ebpmf_identity(X,3,printevery = 1,convergence_criteria = 'mKLabs',ebpm.fn = ebpm::ebpm_point_gamma,smooth_F = F)
plot(fit0$EF[,1],type='l')
plot(fit0$EF[,2],type='l')
plot(fit0$EF[,3],type='l')
fit1 = ebpmf_identity(X,3,printevery = 1,ebpm.fn = c(ebpm::ebpm_point_gamma,smashrgen::ebps),convergence_criteria = 'mKLabs')
plot(fit1$EF_smooth[,1],type='l')
plot(fit1$EF_smooth[,2],type='l')
plot(fit1$EF_smooth[,3],type='l')
