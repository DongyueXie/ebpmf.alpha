#'@title Make Init L and F on similar scale
#'@param L N by K loading matrix
#'@param FF K by p factor matrix
#'@details LF = L \times S_l^{-1} \times S_l \times S_f \times S_f^{-1} \times F
#'
adjLF = function(L,FF){
  gammaL = colSums(L)
  gammaF = rowSums(FF)
  adjScale = sqrt(gammaL*gammaF)
  L = t(t(L) * (adjScale/gammaL))
  FF = FF * (adjScale/gammaF)
  return(list(L_init = L, F_init=FF))
}
