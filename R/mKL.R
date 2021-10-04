#'@title calculate mean KL divergence of 2 nonnegative matrices
#'@export
#'

mKL = function(A,B){
  D = A*log(A/B)-A+B
  mean(as.matrix(D),na.rm=T)
}
