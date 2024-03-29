
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ebpmf.alpha

<!-- badges: start -->
<!-- badges: end -->

An R Package for factorizing count matrix using flexible empirical Bayes
methods.

## Installation

You can install the development version of `ebpmf.alpha` from
[GitHub](https://github.com/) with:

``` r
### install dependencies if necessary
# install.packages("devtools")
# devtools::install_github("DongyueXie/ebpm")
# devtools::install_github("DongyueXie/smashr")
# devtools::install_github("DongyueXie/vebpm")
# devtools::install_github("DongyueXie/smashrgen")

devtools::install_github("DongyueXie/ebpmf.alpha")
```

## Example

### ebpmf_log

``` r
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

fit <- ebpmf_log(Y,l0=0,f0=0,flash_control=list(fix_f0=T))
plot(fit$K_trace)
plot(fitted(fit$fit_flash),tcrossprod(Ltrue,Ftrue),col='grey80')
abline(a=0,b=1)
fit$fit_flash$pve
for(k in 1:fit$fit_flash$n_factors){
  plot(fit$fit_flash$F_pm[,k],type='l')
}
```

### ebpmf_identity

``` r
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

fit1 = ebpmf_identity(X,K)
fit0 = fastTopics::fit_poisson_nmf(X,K)
plot(fit0$F[,1],type='l',ylab='',main='factor 1 unsmoothed')
plot(fit0$F[,2],type='l',ylab='',main='factor 2 unsmoothed')
plot(fit0$F[,3],type='l',ylab='',main='factor 3 unsmoothed')
plot(fit1$EF[,1],type='l',ylab='',main='factor 1 smoothed')
plot(fit1$EF[,2],type='l',ylab='',main='factor 2 smoothed')
plot(fit1$EF[,3],type='l',ylab='',main='factor 3 smoothed')
```
