#'@title Smoothed Poisson Matrix Factorization
#'@param X: count matrix
#'@param K: number of factors/ranks
#'@param init: initialization methods, 'lee','scd' from package NNLM, or 'uniform' randomly initialize; or provide init as a list with L_init and F_init.
#'@param maxiter: maximum iterations
#'@param tol: stop criteria
#'@param bmsm_control: control parameters of BMSM, see bmsm_control_default()
#'@param ebpm_method point_gamma or two_gamma
#'@param ebpm_control: control parameters of ebpm, see ebpm_control_default()
#'@param nug_control control parameters of smashgen, see nug_control_default()
#'@param smooth_f,smooth_l: whether to get smooth estimate of loadings or factors.
#'@param filter.number,family: wavelet basis, see wavethresh pakcage for more details.
#'@param rounding: whether to round Z after each iteration.
#'@return ql,qf: posterior of loadings and factors
#'@return Lambda: estimated lambda; init: initial values
#'@return KL: mean KL of lambda and estimated lambda
#'@import mixsqp
#'@import NNLM
#'@import ebpm
#'@export

stm = function(X,K,init = 'lee',maxiter=100,tol=1e-4,
                  bmsm_control_l=list(), bmsm_control_f=list(),
                  nug_control_l=list(), nug_control_f=list(),
                  filter.number = 1,family = "DaubExPhase",
                  ebpm_method='point_gamma',
                  ebpm_control_l=list(), ebpm_control_f=list(),
                  smooth_f=TRUE,smooth_l=FALSE,rounding=FALSE,
                  nugget=FALSE,printevery=10,return_all=TRUE){

  #initialize q(Z)

  n = dim(X)[1]
  p = dim(X)[2]

  if(is.list(init)){
    L_init = init$L_init
    F_init = init$F_init

  }else{
    if(init == 'scd'){
      X_init_fit = NNLM::nnmf(X,K,method='scd',loss='mkl',show.warning = F,verbose = F,max.iter = 50)
      L_init = X_init_fit$W
      F_init = X_init_fit$H
    }
    if(init == 'lee'){
      X_init_fit = NNLM::nnmf(X,K,method='lee',loss='mkl',show.warning = F,verbose = F,max.iter = 100)
      L_init = X_init_fit$W
      F_init = X_init_fit$H
    }
    if(init == 'uniform'){
      L_init = matrix(runif(n*K),nrow=n,ncol=K)
      F_init = matrix(runif(K*p),nrow=K,ncol=p)
      ratio = mean(X)/(mean(L_init)*mean(F_init))
      L_init = L_init*sqrt(ratio)
      F_init = F_init*sqrt(ratio)
    }
  }

  ratio = adjLF(L_init,F_init)
  L_init = ratio$L_init
  F_init = ratio$F_init


  ql_hat = list(El = L_init, Elogl = log(L_init))
  qf_hat = list(Ef = F_init, Elogf = log(F_init))

  gl_hat = list()
  gf_hat = list()

  EZ = array(dim = c(n,p,K))

  EZ = Calc_EZ(X,K,EZ,ql_hat,qf_hat)

  KL = c()
 # loglik = c()

  KL[1] = mKL(X,L_init%*%F_init)


  for(iter in 1:maxiter){

    if(rounding){
      EZ = round(EZ)
    }

    #loglikL=0
    #loglikR=0

    nugget_l = rep(0,K)
    nugget_f = rep(0,K)

    for(k in 1:K){




      l_seq = rowSums(EZ[,,k])
      l_scale = sum(qf_hat$Ef[k,])


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
        lk_hat = update_smooth(l_seq, l_scale, nugget,bmsm_control_l,nug_control_l,
                               filter.number = filter.number ,family = family)
        ql_hat$El[,k] = lk_hat$E
        ql_hat$Elogl[,k] = lk_hat$Elog
        #loglikL = loglikL + lk_hat$loglik
        nugget_l[k] = lk_hat$nugget

        gl_hat[[k]] = lk_hat$pi_weights

      }else{
        lk_hat = update_nsmooth(l_seq,l_scale,ebpm_control_l,ebpm_method)
        ql_hat$El[,k] = lk_hat$posterior$mean
        ql_hat$Elogl[,k] = lk_hat$posterior$mean_log
        #loglikL = loglikL + lk_hat$log_likelihood

        gl_hat[[k]] = lk_hat$fitted_g
      }


      # Update F

      f_seq = colSums(EZ[,,k])
      f_scale = sum(ql_hat$El[,k])

      if(smooth_f){
        fk_hat = update_smooth(f_seq, f_scale,nugget,bmsm_control_f,nug_control_f,
                               filter.number = filter.number ,family = family)
        qf_hat$Ef[k,] = fk_hat$E
        qf_hat$Elogf[k,] = fk_hat$Elog
        #loglikR = loglikR + fk_hat$loglik
        nugget_f[k] = fk_hat$nugget

        gf_hat[[k]] = fk_hat$pi_weight

      }else{
        fk_hat = update_nsmooth(f_seq,f_scale,ebpm_control_f,ebpm_method)
        qf_hat$Ef[k,] = fk_hat$posterior$mean
        qf_hat$Elogf[k,] = fk_hat$posterior$mean_log
        #loglikR = loglikR + fk_hat$log_likelihood
        gf_hat[[k]] = fk_hat$fitted_g
      }



    }
    # Update Z
    EZ = Calc_EZ(X,K,EZ,ql_hat,qf_hat)

    KL[iter+1] = mKL(X,ql_hat$El%*%qf_hat$Ef)
    #loglik[iter] = loglikR+loglikL

    ########
    if(iter%%printevery==0){
      print(sprintf('At iter %d, mKL: %f',iter,KL[iter+1]))
    }
    ########

    if(abs(KL[iter+1]-KL[iter])<=tol){
      break
    }

  }
  if(iter==maxiter){
    warning('Reached maximum iterations')
  }

  lambda_hat = ql_hat$El%*%qf_hat$Ef
  lambda_init = L_init%*%F_init

  loglik = sum(dpois(X,lambda_hat,log = TRUE))

  if(return_all){
    return(list(ql=ql_hat,qf=qf_hat,gf=gf_hat,gl=gl_hat,KL=KL,loglik=loglik,Lambda_hat=lambda_hat,
                init = list(L_init = L_init,F_init=F_init,Lambda_init = lambda_init),EZ=EZ,
                input = list(X=X,K=K),nugget=list(nugget_l=nugget_l,nugget_f=nugget_f)))
  }else{
    return(list(ql=ql_hat$El,qf=qf_hat$Ef,nugget=list(nugget_l=nugget_l,nugget_f=nugget_f),KL=KL))
  }

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

update_smooth = function(x,sf,nugget,bmsm_control=list(),nug_control=list(),
                         filter.number = 1,family = "DaubExPhase"){
  if(min(x) < 0){stop ("negative values in x not permitted")}
  if(nugget){

    control0 = nug_control_default()
    if (any(!is.element(names(nug_control),names(control0))))
      stop("Argument \"nug_control\" contains unknown parameter names")

    control1 = modifyList(control0,nug_control,keep.null = TRUE)

    fit = smash.gen.poiss(x,s=sf,filter.number=filter.number,
                          family=family,nugget=control1$nugget,
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
    out = ebpm_point_gamma(x,s,g_init,fix_g,pi0,control)
  }
  if(ebpm_method=='two_gamma'){
    out = ebpm_two_gamma(x,s,g_init,fix_g,pi0,control)
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
       tol=1e-2)
}


