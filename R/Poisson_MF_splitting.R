#'@title Fit Sparse Poisson matrix factorization
#'@param Y count data matrix
#'@param S The known scaling factor matrix, background frequency.
#'@param sigma2 the variance term
#'@param est_sigma2 whether estimate the variance term or fix it
#'@return fitted object
#'@import flashier
#'@importFrom vebpm pois_mean_GG
#'@export
splitting_PMF_flashier = function(Y,S,sigma2=NULL,est_sigma2 = TRUE,
                         Kmax=10,var_type='by_row',
                         M_init = NULL,
                         maxiter=100,tol=0.1,
                         verbose_flash=0,
                         printevery=10){

  n = nrow(Y)
  p = ncol(Y)

  if(is.null(S)){
    S = 1
  }
  if(length(S)==1){
    S = matrix(S,nrow=n,ncol=p)
  }

  if(is.null(sigma2)|is.null(M_init)){
    # pre-estimate sigma2, assuming LF = 0?.
    if(var_type=='constant'){
      init_val = pois_mean_GG(as.vector(Y),as.vector(S),prior_mean = 0,prior_var = NULL,tol=1e-3)
      sigma2_init = init_val$fitted_g$var
      M0 = matrix(init_val$posterior$mean_log,nrow=n,ncol=p)
    }
    if(var_type=='by_row'){
      init_val = lapply(1:n,function(i){
        return(pois_mean_GG(Y[i,],S[i,],prior_mean = 0,prior_var = NULL,tol=1e-3))
      })
      sigma2_init = unlist(lapply(init_val,function(fit){fit$fitted_g$var}))
      M0 = do.call(rbind,lapply(init_val,function(fit){fit$posterior$mean_log}))
    }
    if(var_type=='by_col'){
      init_val = lapply(1:p,function(i){
        return(pois_mean_GG(Y[,i],S[,i],prior_mean = 0,prior_var = NULL,tol=1e-3))
      })
      sigma2_init = unlist(lapply(init_val,function(fit){fit$fitted_g$var}))
      M0 = do.call(cbind,lapply(init_val,function(fit){fit$posterior$mean_log}))
    }
  }
  if(is.null(sigma2)){
    sigma2 = sigma2_init
    est_sigma2 = TRUE
    rm(sigma2_init)
  }
  if(is.null(M_init)){
    M = M0
    rm(M0)
    rm(init_val)
  }else{
    M = M_init
  }

  const = sum(Y*log(S)) - sum(lfactorial(Y))
  ## fit flashier on M, sigma2 with backfitting
  if(var_type=='by_row'){
    S.dim = 1
  }else if(var_type=='by_col'){
    S.dim = 2
  }else if(var_type=='constant'){
    S.dim = NULL
  }else{
    stop('Non-supported var type')
  }

  fit_flash = flash.init(M, S = sqrt(sigma2), var.type = NULL, S.dim = S.dim)%>%
    flash.add.greedy(Kmax = Kmax,verbose = verbose_flash) %>%
    flash.backfit(verbose = verbose_flash) %>%
    flash.nullcheck(verbose = verbose_flash)

  KL_LF = sum(ff.KL(fit_flash$flash.fit,1)) + sum(ff.KL(fit_flash$flash.fit,2))
  V = matrix(1/n,nrow=n,ncol=p)
  obj = calc_split_PMF_obj_flashier(Y,S,sigma2,M,V,fit_flash,KL_LF,const,var_type)

  start_time = Sys.time()

  for(iter in 1:maxiter){

    res = vga_pois_solver_mat(M,Y,S,fitted(fit_flash),adjust_var_shape(sigma2,var_type,n,p))
    M = res$M
    V = res$V
    #print(paste('After vga,elbo is',calc_split_PMF_obj_flashier(Y,S,sigma2,M,V,fit_flash,KL_LF,const,var_type)))

    # update sigma2
    if(est_sigma2){
      if(var_type=='constant'){
        sigma2 = (sum(V) +sum(fit_flash$flash.fit$R2))/(n*p)
      }else if(var_type=='by_row'){
        sigma2 = (rowSums(V)+fit_flash$flash.fit$R2)/p
      }else if(var_type=='by_col'){
        sigma2 = (colSums(V)+fit_flash$flash.fit$R2)/n
      }else{
        stop('Non-supported var type')
      }
    }
    #print(paste('After sigma2,elbo is',calc_split_PMF_obj_flashier(Y,S,sigma2,M,V,fit_flash,KL_LF,const,var_type)))
    ## solve flash
    fit_flash = flash.init(M, S = sqrt(sigma2), var.type = NULL, S.dim = S.dim)%>%
      flash.init.factors(init = fit_flash) %>%
      flash.add.greedy(Kmax = Kmax,verbose = verbose_flash) %>%
      flash.backfit(verbose = verbose_flash) %>%
      flash.nullcheck(verbose = verbose_flash)

    # fit_flash = flash.init(M, S = sqrt(sigma2), var.type = NULL, S.dim = S.dim)%>%
    #   flash.add.greedy(Kmax = Kmax,verbose = verbose_flash)


    KL_LF = sum(ff.KL(fit_flash$flash.fit,1)) + sum(ff.KL(fit_flash$flash.fit,2))
    #print(paste('After flash,elbo is',calc_split_PMF_obj_flashier(Y,S,sigma2,M,V,fit_flash,KL_LF,const,var_type)))




    # check convergence
    obj[iter + 1] = calc_split_PMF_obj_flashier(Y,S,sigma2,M,V,fit_flash,KL_LF,const,var_type)
    if((obj[iter+1] - obj[iter]) < tol){
      if((obj[iter+1] - obj[iter])<0){
        warning('An iteration decreases ELBO')
      }
      break
    }

    if(iter%%printevery==0){
      print(paste('At iter',iter, 'ELBO=',obj[iter+1]))
    }

  }

  end_time = Sys.time()
  return(list(fit_flash=fit_flash,elbo=obj[length(obj)],eblo_trace=obj,
              sigma2 = sigma2,run_time = difftime(end_time,start_time,units='auto'),M=M,V=V))
}



calc_split_PMF_obj_flashier = function(Y,S,sigma2,M,V,fit_flash,KL_LF,const,var_type){
  R2 = fit_flash$flash.fit$R2
  n = nrow(Y)
  p = ncol(Y)
  if(var_type=='by_row'){
    sv = rowSums(V)
    ss = p
  }else if(var_type=='by_col'){
    sv = colSums(V)
    ss = n
  }else if(var_type=='constant'){
    sv = sum(V)
    ss = n*p
  }else{
    stop('Non-supported var type')
  }
  val = sum(Y*M - S*exp(M+V/2)   + log(2*pi*V)/2 + 0.5 ) - sum(ss*log(2*pi*sigma2)/2)- sum(sv/2/sigma2) - sum(R2/2/sigma2) + const+ KL_LF
  return(val)
}


adjust_var_shape = function(sigma2,var_type,n,p){
  if(var_type=='constant'){
    sigma2 = matrix(sigma2,nrow=n,ncol=p)
  }else if(var_type=='by_row'){
    sigma2 = matrix(sigma2,nrow=n,ncol=p,byrow = F)
  }else if(var_type=='by_col'){
    sigma2 = matrix(sigma2,nrow=n,ncol=p,byrow = T)
  }else{
    stop('Non-supported var type')
  }
  sigma2
}


#'@title a matrix version of the solver
#'@importFrom vebpm vga_pois_solver
vga_pois_solver_mat = function(init_Val,X,S,Beta,Sigma2,maxiter=1000,tol=1e-8){

  n = nrow(X)
  p = ncol(X)
  # transform them to vector
  x = as.vector(X)
  s = as.vector(S)
  beta = as.vector(Beta)
  sigma2 = as.vector(Sigma2)
  init_val = as.vector(init_Val)

  res = vga_pois_solver(init_val,x,s,beta,sigma2,maxiter=maxiter,tol=tol,method='newton')

  return(list(M = matrix(res$m,nrow=n,ncol=p,byrow = F),V = matrix(res$v,nrow=n,ncol=p,byrow = F)))

}








