sla_full <- readRDS("~/Downloads/sla_full.rds")
dim(sla_full$data)
sum(sla_full$data==0)/prod(dim(sla_full$data))

doc_to_use = order(rowSums(sla_full$data),decreasing = T)[1:round(nrow(sla_full$data)*0.6)]
mat = sla_full$data[doc_to_use,]
samples = sla_full$samples
samples = lapply(samples, function(z){z[doc_to_use]})
word_to_use = which(colSums(mat>0)>=5)
mat = mat[,word_to_use]

set.seed(1)
library(ebpmf)
fit_ebpmf_K1 = ebpmf_log(mat,
                         flash_control=list(backfit_extrapolate=T,backfit_warmstart=T,
                                            ebnm.fn = c(ebnm::ebnm_point_exponential, ebnm::ebnm_point_exponential),
                                            loadings_sign = 1,factors_sign=1,Kmax=1),
                         init_control = list(n_cores=5,flash_est_sigma2=F,log_init_for_non0y=T),
                         general_control = list(maxiter=100,save_init_val=T,save_latent_M=T),
                         sigma2_control = list(return_sigma2_trace=T))
fit_ebpmf_K1$fit_flash$L_pm
plot(fit_ebpmf_K1$fit_flash$F_pm[,2])
resid  = flashier:::residuals.flash(fit_ebpmf_K1$fit_flash)

fit_ebpmf_Kmax = ebpmf_log(mat,
                         flash_control=list(backfit_extrapolate=T,backfit_warmstart=T,
                                            ebnm.fn = c(ebnm::ebnm_point_exponential, ebnm::ebnm_point_exponential),
                                            loadings_sign = 1,factors_sign=1,Kmax=100),
                         init_control = list(n_cores=5,flash_est_sigma2=T,log_init_for_non0y=T),
                         general_control = list(maxiter=100,save_init_val=T,save_latent_M=T),
                         sigma2_control = list(return_sigma2_trace=T))
L= fit_ebpmf_Kmax$fit_flash$L_pm[,-c(1,2)]
F_pm = fit_ebpmf_Kmax$fit_flash$F_pm[,-c(1,2)]
rownames(L)<-1:nrow(L)

Lnorm = t(t(L)/apply(L,2,max))
Fnorm = t(t(F_pm)*apply(L,2,max))
khat = apply(Lnorm,1,which.max)
Lmax = apply(Lnorm,1,max)
plot(Lmax)

khat[Lmax<0.1] = 0
keyw.nn =list()

for(k in 1:ncol(Fnorm)){
  key = Fnorm[,k]>log(2)
  keyw.nn[[k]] = (colnames(mat)[key])[order(Fnorm[key,k],decreasing = T)]
}
print(keyw.nn)

structure_plot_general = function(Lhat,Fhat,grouping,title=NULL,
                                  loadings_order = 'embed',
                                  print_plot=FALSE,
                                  seed=12345,
                                  n_samples = NULL,
                                  gap=40,
                                  std_L_method = 'sum_to_1',
                                  show_legend=TRUE,
                                  K = NULL
){
  set.seed(seed)
  #s       <- apply(Lhat,2,max)
  #Lhat    <-   t(t(Lhat) / s)

  if(is.null(n_samples)&all(loadings_order == "embed")){
    n_samples = 2000
  }

  if(std_L_method=='sum_to_1'){
    Lhat = Lhat/rowSums(Lhat)
  }
  if(std_L_method=='row_max_1'){
    Lhat = Lhat/c(apply(Lhat,1,max))
  }
  if(std_L_method=='col_max_1'){
    Lhat = apply(Lhat,2,function(z){z/max(z)})
  }
  if(std_L_method=='col_norm_1'){
    Lhat = apply(Lhat,2,function(z){z/norm(z,'2')})
  }

  if(!is.null(K)){
    Lhat = Lhat[,1:K]
    Fhat = Fhat[,1:K]
  }
  Fhat = matrix(1,nrow=3,ncol=ncol(Lhat))
  if(is.null(colnames(Lhat))){
    colnames(Lhat) <- paste0("k",1:ncol(Lhat))
  }
  fit_list     <- list(L = Lhat,F = Fhat)
  class(fit_list) <- c("multinom_topic_model_fit", "list")
  p <- structure_plot(fit_list,grouping = grouping,
                      loadings_order = loadings_order,
                      n = n_samples,gap = gap,verbose=F) +
    labs(y = "loading",color = "dim",fill = "dim") + ggtitle(title)
  if(!show_legend){
    p <- p + theme(legend.position="none")
  }
  if(print_plot){
    print(p)
  }
  return(p)
}

structure_plot_general(Lnorm,Fnorm,grouping = samples$journal,std_L_method = 'col_max_1')
