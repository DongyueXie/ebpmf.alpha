
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ebpmf

<!-- badges: start -->
<!-- badges: end -->

An R Package.

## Installation

You can install the development version of `ebpmf` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
# devtools::install_github("DongyueXie/vebpm")
# devtools::install_github("DongyueXie/smashrgen")
# devtools::install_github("DongyueXie/smashr")
# devtools::install_github("DongyueXie/ebpm")
devtools::install_github("DongyueXie/ebpmf")
```

## Example

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
