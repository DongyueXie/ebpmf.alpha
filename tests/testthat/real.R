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

S = tcrossprod(c(rowSums(Y)),c(colSums(Y)))/sum(Y)
Y = as.matrix(Y)
# run method
fit = splitting_PMF_flashier(Y,S,var_type = 'by_col',Kmax = 30,maxiter = 1000,verbose = TRUE)

