datax = read.csv('/project2/mstephens/gtex-stm/Counts/TPM3.Counts.csv.gz')
datax = datax[,-1]
datax = as.matrix(datax)
library(stm)
fit = ebpmf_identity(datax,K=7,printevery = 1)


datax = read.csv('/project2/mstephens/gtex-stm/Counts/SRSF3.Counts.csv.gz')
datax = datax[,-1]
datax = as.matrix(datax)
library(stm)
fit = ebpmf_identity(datax,K=7,printevery = 1)
