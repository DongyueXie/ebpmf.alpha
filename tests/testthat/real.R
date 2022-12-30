library(fastTopics)
library(Matrix)
data(pbmc_facs)
counts <- pbmc_facs$counts
table(pbmc_facs$samples$subpop)
## use only B cell and NK cell
cells = pbmc_facs$samples$subpop%in%c('B cell', 'NK cell')
Y = counts[cells,]
dim(Y)
# filter out genes that has few expressions(3% cells)
genes = (colSums(Y>0) > 0.03*dim(Y)[1])
Y = Y[,genes]
# make sure there is no zero col and row
sum(rowSums(Y)==0)
sum(colSums(Y)==0)
dim(Y)

rm(counts)
rm(genes)
rm(cells)

#Y = as.matrix(pbmc_facs$counts)
S = tcrossprod(c(rowSums(Y)),c(colSums(Y)))/sum(Y)
Y = as.matrix(Y)
# run method
fit = splitting_PMF_flashier(Y,S,var_type = 'by_col',Kmax = 30,maxiter = 1000,verbose = TRUE,n_cores = 10)

##############
fit = readRDS('/project2/mstephens/dongyue/poisson_mf/pbmc_3cells_Sij.rds')
kset = order(fit$fit$fit_flash$pve,decreasing = TRUE)[1:3]
Sigma2 = matrix(fit$fit$sigma2,nrow=nrow(fit$S),ncol=ncol(fit$S),byrow = T)
simdata = sim_data_real(fit$S,fit$fit$fit_flash$L.pm[,kset],fit$fit$fit_flash$F.pm[,kset],Sigma2,n_simu = 2,seed=12345)
rm(Sigma2)
rm(fit)
out = simu_study_PMF(simdata,n_cores = 1,Kmax=30,var_type='by_col',maxiter=100,tol=1e-5)
saveRDS(out,file=paste('PMF',n_simu,'_n',n,'_K',K,'_pbmc_3cells','.rds',sep=''))
