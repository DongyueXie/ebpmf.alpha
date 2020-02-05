#'@title calculate mean KL divergence of 2 nonnegative matrices
#'@export
#'

mKL = function(A,B){
  D = A*log(A/B)-A+B
  D[is.nan(D)] = 0
  mean(D)
}
