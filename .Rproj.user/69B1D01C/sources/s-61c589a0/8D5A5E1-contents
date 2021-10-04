library(stm)
library(ebpm)
library(NNLM)
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

lambda = L %*% FF
X = matrix(rpois(n=length(lambda),lambda),nrow=n)
image(X)
X = Matrix::Matrix(X,sparse = T)
fit_smooth = stm(X,3,nugget = F,tol=1e-4)

plot(fit_smooth$EF[,1],type='l')
plot(fit_smooth$EF[,2],type='l')
plot(fit_smooth$EF[,3],type='l')

fit_smooth = stm(X,3,nugget = T,tol=1e-4)

plot(fit_smooth$EF[,1],type='l')
plot(fit_smooth$EF[,2],type='l')
plot(fit_smooth$EF[,3],type='l')


