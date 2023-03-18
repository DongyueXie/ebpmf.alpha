datax = read.csv('/project2/mstephens/gtex-stm/Counts/TPM3.Counts.csv.gz')
datax = datax[,-1]
datax = as.matrix(datax)
fit = ebpmf_identity(datax,K=7,printevery = 1)


datax = read.csv('/project2/mstephens/gtex-stm/Counts/SRSF3.Counts.csv.gz')
datax = datax[,-1]
datax = as.matrix(datax)
fit = ebpmf_identity(datax,K=7,printevery = 1)

tpm3 = tpm3[,colSums(tpm3)>0]
fit = ebpmf_identity(tpm3,7,init = list(L_init = tpm3_fit_fasttopics$Ln,F_init = tpm3_fit_fasttopics$Fn),printevery = 1)
