library(flashier)
data("gtex")
fl <- flash.init(gtex,S=rep(1,ncol(gtex)),var.type = NULL) %>%
  flash.add.greedy(Kmax = 3) %>%
  flash.backfit(maxiter=1) %>%
  flash.nullcheck()


temp = flash.add.greedy(fl,Kmax=1,verbose = 0)

temp$flash.fit$Y = gtex+rnorm(prod(dim(gtex)))
temp$residuals.sd = sqrt(rep(2,ncol(gtex)))
temp$flash.fit$given.tau = rep(1/2,ncol(gtex))
temp$flash.fit$est.tau = rep(1/2,ncol(gtex))
temp$flash.fit$tau = rep(1/2,ncol(gtex))

temp2 = flash.add.greedy(temp,Kmax=1,verbose = 2)


ones <- matrix(1, nrow = nrow(gtex), ncol = 1)
ls.soln <- t(solve(crossprod(ones), crossprod(ones, gtex)))
flinit <- flash.init(gtex) %>%
  flash.init.factors(init = list(ones, ls.soln)) %>%
  flash.add.greedy(Kmax=1) %>%
  flash.backfit(maxiter = 1)

fit = flashier:::set.flash.data(gtex, S = rep(1,ncol(gtex)), S.dim = 2, var.type = NULL)
