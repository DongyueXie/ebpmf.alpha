#'@title standard loadings and factors from Poisson matrix Factorization
#'@param L: n by k matrix
#'@param FF: k by p matrix
#'@export
#'

poisson_to_multinom <- function (FF, L) {
  L <- t(t(L) * colSums(FF))
  s <- rowSums(L)
  L <- L / s
  FF <- scale.cols(FF)
  return(list(FF = FF,L = L,s = s))
}

scale.cols <- function (A)
  apply(A,2,function (x) x/sum(x))
