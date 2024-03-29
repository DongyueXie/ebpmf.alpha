#'@title Smoothed Poisson Matrix Factorization
#'@description Smoothed Poisson Matrix Factorization/ smoothed topic model
#'@param X count matrix
#'@param K number of factors/ranks
#'@param init initialization methods, 'lee','scd' from package NNLM, or 'uniform' randomly initialize; or provide init as a list with L_init and F_init.
#'@param init_loss loss function of the initialization method, either mkl or mse.
#'@param maxiter maximum iterations
#'@param tol stop criteria
#'@param fix_F if TRUE, F will no be updated.
#'@param bmsm_control control parameters of BMSM, see bmsm_control_default()
#'@param ebpm_method point_gamma or two_gamma
#'@param ebpm_control control parameters of ebpm, see ebpm_control_default()
#'@param nug_control control parameters of smashgen, see nug_control_default()
#'@param smooth_f,smooth_l whether to get smooth estimate of loadings or factors.
#'@param nugget whether to assume nugget effects
#'@param return_all whether return all outputs or simplified ones
#'@return EL,EF: posterior of loadings and factors
#'@import mixsqp
#'@import NNLM
#'@import ebpm
#'@import Matrix
#'@export

stm = function(X,K,
               init = 'scd',init_loss = 'mkl',maxiter=100,tol=1e-3,
               fix_F = FALSE,
                  bmsm_control_l=list(), bmsm_control_f=list(),
                  nug_control_l=list(), nug_control_f=list(),
                  #filter.number = 1,family = "DaubExPhase",
                  ebpm_method='point_gamma',
                  ebpm_control_l=list(), ebpm_control_f=list(),
                  smooth_f=TRUE,smooth_l=FALSE,
                  nugget=FALSE,
               printevery=10){


  n = dim(X)[1]
  p = dim(X)[2]

  res = init_stm(X,K,init,init_loss)
  #plot(res$ql$El[,1])
  #plot(res$ql$El[,2])
  #plot(res$ql$El[,3])
  #inited = list(L_init = res$ql$El,F_init = res$qf$Ef)

  #EZ = array(dim = c(n,p,K))

  #EZ = Calc_EZ(X,K,EZ,ql_hat,qf_hat)

  KL = c()
 # loglik = c()

  #browser()
  KL[1] = mKL(X,tcrossprod(res$ql$El,res$qf$Ef))

  X = Matrix::Matrix(X,sparse = TRUE)
  X_idx = summary(X)

  for(iter in 1:maxiter){

    b_k_max = 0

    for(k in 1:K){

      # get row and col sums of EZ_k
      b_k = res$ql$Elogl[X_idx$i,k]+res$qf$Elogf[X_idx$j,k] - res$a

      EZ_k = sparseMatrix(i=X_idx$i,j=X_idx$j,x = X_idx$x*exp(b_k)/res$b,dims = c(n,p))


      l_seq = rowSums(EZ_k)
      l_scale = sum(res$qf$Ef[,k])


      # adj.ratio = sqrt(l_scale/f_scale)
      #
      # l_scale = l_scale/adj.ratio
      # l_seq = l_seq * adj.ratio
      # f_scale = f_scale*adj.ratio
      # f_seq = f_seq / adj.ratio


      #print(l_scale)
      #print(f_scale)

      # Update L
      if(smooth_l){
        lk_hat = update_smooth(l_seq, l_scale, nugget,bmsm_control_l,nug_control_l)
        res$ql$El[,k] = lk_hat$E
        res$ql$Elogl[,k] = lk_hat$Elog
        #loglikL = loglikL + lk_hat$loglik
        res$nugget_l[k] = lk_hat$nugget

        res$gl[[k]] = lk_hat$pi_weights

      }else{
        lk_hat = update_nsmooth(l_seq,l_scale,ebpm_control_l,ebpm_method)
        res$ql$El[,k] = lk_hat$posterior$mean
        res$ql$Elogl[,k] = lk_hat$posterior$mean_log
        #loglikL = loglikL + lk_hat$log_likelihood
        res$gl[[k]] = lk_hat$fitted_g
      }


      # Update F

      if(!fix_F){
        f_seq = colSums(EZ_k)
        f_scale = sum(res$ql$El[,k])

        if(smooth_f){
          fk_hat = update_smooth(f_seq, f_scale,nugget,bmsm_control_f,nug_control_f)
          res$qf$Ef[,k] = fk_hat$E
          res$qf$Elogf[,k] = fk_hat$Elog
          #loglikR = loglikR + fk_hat$loglik
          res$nugget_f[k] = fk_hat$nugget
          res$gf[[k]] = fk_hat$pi_weight

        }else{
          fk_hat = update_nsmooth(f_seq,f_scale,ebpm_control_f,ebpm_method)
          res$qf$Ef[,k] = fk_hat$posterior$mean
          res$qf$Elogf[,k] = fk_hat$posterior$mean_log
          #loglikR = loglikR + fk_hat$log_likelihood
          res$gf[[k]] = fk_hat$fitted_g
        }
      }

      b_k_new = res$ql$Elogl[X_idx$i,k] + res$qf$Elogf[X_idx$j, k] - res$a
      res$b = res$b - exp(b_k) + exp(b_k_new)
      b_k_max = pmax(b_k_new, b_k_max)


    }
    # # Update Z
    # EZ = Calc_EZ(X,K,EZ,ql_hat,qf_hat)

    res$b = res$b/exp(b_k_max)
    res$a = b_k_max + res$a


    KL[iter+1] = mKL(X,tcrossprod(res$ql$El,res$qf$Ef))
    #loglik[iter] = loglikR+loglikL

    ########
    if(iter%%printevery==0){
      print(sprintf('At iter %d, mean KL: %f',iter,KL[iter+1]))
    }
    ########

    if(abs(KL[iter+1]-KL[iter])<=tol){
      break
    }

  }
  if(iter==maxiter){
    warning('Reached maximum iterations')
  }

  lambda_hat = tcrossprod(res$ql$El,res$qf$Ef)
  #lambda_init = L_init%*%F_init

  # loglik = sum(dpois(X,lambda_hat,log = TRUE))

  ldf = poisson2multinom(res$qf$Ef,res$ql$El)
  fit = list(res = res,EL = ldf$L,EF = ldf$FF,d=ldf$s)
  return(fit)
  # if(return_all){
  #   return(list(ql=ql_hat,qf=qf_hat,gf=gf_hat,gl=gl_hat,KL=KL,Lambda_hat=lambda_hat,
  #               init = inited,
  #               input = list(X=X,K=K),nugget=list(nugget_l=nugget_l,nugget_f=nugget_f)))
  # }else{
  #   return(list(ql=ql_hat$El,qf=qf_hat$Ef,nugget=list(nugget_l=nugget_l,nugget_f=nugget_f),KL=KL))
  # }

}

#'@title initialize the stm model
#'@param X input data matrix
#'@param K number of topics
#'@param init init methods, or a list of init L and F
#'@param init_loss mkl or mse
#'@export
init_stm = function(X,K,init,init_loss){

  if(is.list(init)){
    L_init = init$L_init
    F_init = init$F_init

    if(is.null(L_init)){
      X_init_fit = NNLM::nnmf(as.matrix(X),K,method='lee',
                              loss='mse',show.warning = F,
                              init = list(H=t(F_init)),
                              verbose = F,max.iter = 50)
      L_init = X_init_fit$W
    }

  }else{

    if(init%in%c('scd','lee')){
      X_init_fit = NNLM::nnmf(as.matrix(X),K,method=init,loss=init_loss,show.warning = F,verbose = F,max.iter = 50)
      L_init = X_init_fit$W
      F_init = t(X_init_fit$H)
    }

    # if(init == 'scd'){
    #   X_init_fit = NNLM::nnmf(as.matrix(X),K,method='scd',loss=init_loss,show.warning = F,verbose = F,max.iter = 50)
    #   L_init = X_init_fit$W
    #   F_init = t(X_init_fit$H)
    # }
    # if(init == 'lee'){
    #   X_init_fit = NNLM::nnmf(as.matrix(X),K,method='lee',loss=init_loss,show.warning = F,verbose = F,max.iter = 100)
    #   L_init = X_init_fit$W
    #   F_init = t(X_init_fit$H)
    # }
    if(init == 'uniform'){
      L_init = matrix(runif(n*K),nrow=n,ncol=K)
      F_init = matrix(runif(K*p),nrow=p,ncol=K)
      ratio = median(X)/(median(L_init)*median(F_init))
      L_init = L_init*sqrt(ratio)
      F_init = F_init*sqrt(ratio)
    }
    if(init == 'kmeans'){
      kmeans.init=kmeans(as.matrix(X),K,nstart=5)
      L_init = rep(1,n)%o%normalize(as.vector(table(kmeans.init$cluster)))
      F_init = t(kmeans.init$centers)
      row.names(F_init)=NULL
    }
  }

  # adjust scale of L and F, mainly for stability.
  ratio = adjLF(L_init,F_init)
  L_init = ratio$L_init
  F_init = ratio$F_init

  Elogl = log(L_init+1e-10)
  Elogf = log(F_init+1e-10)

  ql = list(El = L_init, Elogl = Elogl)
  qf = list(Ef = F_init, Elogf = Elogf)

  a = 0
  b = 0

  X = Matrix(X,sparse = TRUE)
  d = summary(X)

  temp = Elogl[d$i,] + Elogf[d$j,]
  a = rowMax(temp)
  b = rowSums(exp(temp-a))


  gl = list()
  gf = list()

  return(list(ql=ql,qf=qf,gl=gl,gf=gf,
              a=a,b=b,
              nugget_l = rep(0,K),
              nugget_f = rep(0,K)))

}

rowMax = function(X){
  do.call(pmax.int, c(na.rm = TRUE, as.data.frame(X)))
}

Calc_EZ = function(X,K,EZ,ql_hat,qf_hat){
  n = nrow(X)
  p = ncol(X)
  for(k in 1:K){
    EZ[,,k] = outer(ql_hat$Elogl[,k], qf_hat$Elogf[k,], "+")
  }
  EZ = softmax3d(EZ)
  EZ = as.vector(EZ)*as.vector(X)
  dim(EZ) = c(n,p,K)
  EZ
}

softmax3d <- function(x){
  score.exp <- exp(x)
  probs <-as.vector(score.exp)/as.vector(rowSums(score.exp,dims=2))
  probs[is.na(probs)] = 0
  dim(probs) <- dim(x)
  return(probs)
}

update_smooth = function(x,sf,nugget,bmsm_control=list(),nug_control=list()){
  if(min(x) < 0){stop ("negative values in x not permitted")}
  if(nugget){

    control0 = nug_control_default()
    if (any(!is.element(names(nug_control),names(control0))))
      stop("Argument \"nug_control\" contains unknown parameter names")

    control1 = modifyList(control0,nug_control,keep.null = TRUE)

    fit = smash.gen.poiss(x,s=sf,filter.number=control1$filter.number,
                          family=control1$family,nugget=control1$nugget,
                          robust=control1$robust,
                          robust.q = control1$robust.q,
                          transformation = control1$transformation,
                          method = control1$method,
                          nug.init = control1$nug.init,
                          ash.pm = control1$ash.pm,
                          eps = control1$eps,
                          maxiter = control1$maxiter,
                          tol = control1$tol)
    est = fit$lambda.est
    if(control1$transformation=='lik_expansion'){
      est_log = fit$mu.est
    }else{
      est_log = log(est)
    }
    pi_weights = NULL
    nugget.est = fit$nugget.est
  }else{
    fit = BMSM(x,sf,bmsm_control)
    est = fit$E
    est_log = fit$Elog
    pi_weights = fit$pi_weights
    nugget.est = 0
  }

  #loglik = fit$loglik

  results = list("E" = est,
                 'Elog' = est_log,
                 "pi_weights" = pi_weights,
                 #"loglik" = loglik,
                 "nugget" = nugget.est)

  return(results)

}


update_nsmooth = function(x,s,ebpm_control = list(),ebpm_method){

  control0 = ebpm_control_default()
  if (any(!is.element(names(ebpm_control),names(control0))))
    stop("Argument \"ebpm_control\" contains unknown parameter names")

  control1 = modifyList(control0,ebpm_control,keep.null = TRUE)

  #scale = control1$scale
  #point_mass=control1$point_mass
  #nullweight=control1$nullweight
  #shape= control1$shape
  g_init = control1$g_init
  fix_g = control1$fix_g
  #m = control1$m
  control =  control1$control
  #low = control1$low
  #d = control1$d
  pi0 = control1$pi0

  if(ebpm_method=='point_gamma'){
    out = ebpm::ebpm_point_gamma(x,s,g_init,fix_g,pi0,control)
  }
  if(ebpm_method=='two_gamma'){
    out = ebpm::ebpm_two_gamma(x,s,g_init,fix_g,pi0,control)
  }

  out



}

#'@title Default parameters of ebpm
#'@export
ebpm_control_default = function(){
  list(pi0 = 'estimate',
       #point_mass=F,
       #nullweight=100,
       #shape=1,
       g_init = NULL,
       fix_g = FALSE,
       #m = 2,
       control =  NULL)
       #low = NULL,
       #d=NULL
}


#'@title Default parameters of smash gen
#'@param filter.number,family wavelet basis, see wavethresh pakcage for more details.
#'@export
nug_control_default = function(){
  list(nugget = NULL,
       robust=T,
       robust.q = 0.99,
       transformation = 'lik_expansion',
       method='ti.thresh',
       nug.init = NULL,
       ash.pm=FALSE,
       eps='estimate',
       maxiter=10,
       tol=1e-2,
       filter.number = 1,
       family = "DaubExPhase")
}


