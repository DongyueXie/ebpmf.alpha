library(stm)
library(NNLM)
library(Matrix)

files = list.files('/project2/mstephens/dongyue/luis/luis')
genes = c()
for(i in 1:length(files)){
  genes[i] = strsplit(files[i],split = '_')[[1]][1]
}
genes = unique(genes)
genes

library(readr)
merge_len = 10

for(g in 8:length(genes)){
  gene = genes[g]
  print(paste('running ',gene))
  RNAseq = read_csv(paste("/project2/mstephens/dongyue/luis/luis/",gene,"_RNAseq.csv.gz",sep=''))
  ATACseq = read_csv(paste("/project2/mstephens/dongyue/luis/luis/",gene,"_ATACseq.csv.gz",sep=''))
  H3K4seq = read_csv(paste("/project2/mstephens/dongyue/luis/luis/",gene,"_H3K4me1.csv.gz",sep=''))

  indis_ATAC = (colnames(ATACseq))[-c(1,2)]
  indis_RNA = (colnames(RNAseq))[-c(1,2)]
  indis_H3K4me1 = (colnames(H3K4seq))[-c(1,2)]
  for(i in 1:length(indis_H3K4me1)){
    name_i = strsplit(indis_H3K4me1[i],split = '_')[[1]]
    indis_H3K4me1[i] = paste(name_i[1],name_i[2],sep = '_')
  }
  indis = intersect(intersect(indis_ATAC,indis_RNA),indis_H3K4me1)

  Y_RNA = t(RNAseq[,match(indis,indis_RNA)+2])
  Y_H3K4 = t(H3K4seq[,match(indis,indis_H3K4me1)+2])
  Y_ATAC = t(ATACseq[,match(indis,indis_ATAC)+2])

  pdf(paste('output/luis/',gene,'_RNAseq.pdf',sep=''),width = 13,height = 5)
  plot(colSums(Y_RNA),type='h',main=paste(gene,'RNAseq'),xlab='base',ylab='counts')
  dev.off()
  pdf(paste('output/luis/',gene,'_H3K4seq.pdf',sep=''),width = 13,height = 5)
  plot(colSums(Y_H3K4),type='h',main=paste(gene,'H3K4me1'),xlab='base',ylab='counts')
  dev.off()
  pdf(paste('output/luis/',gene,'_ATACseq.pdf',sep=''),width = 13,height = 5)
  plot(colSums(Y_ATAC),type='h',main=paste(gene,'ATACseq'),xlab='base',ylab='counts')
  dev.off()

  if(var(c(ncol(Y_RNA),ncol(Y_H3K4),ncol(Y_ATAC)))!=0){
    l = min(c(ncol(Y_RNA),ncol(Y_H3K4),ncol(Y_ATAC)))
    Y_RNA = Y_RNA[,1:l]
    Y_H3K4 = Y_H3K4[,1:l]
    Y_ATAC = Y_ATAC[,1:l]
  }

  if(ncol(Y_RNA)%%merge_len!=0){
    l = ncol(Y_RNA)
    l = floor(l/merge_len)*merge_len
    Y_RNA = Y_RNA[,1:l]
    Y_H3K4 = Y_H3K4[,1:l]
    Y_ATAC = Y_ATAC[,1:l]
  }

  Y_RNAr = do.call(cbind, by(t(Y_RNA), (seq(ncol(Y_RNA)) - 1) %/% merge_len, FUN = colSums))
  Y_H3K4r = do.call(cbind, by(t(Y_H3K4), (seq(ncol(Y_H3K4)) - 1) %/% merge_len, FUN = colSums))
  Y_ATACr = do.call(cbind, by(t(Y_ATAC), (seq(ncol(Y_ATAC)) - 1) %/% merge_len, FUN = colSums))

  X = cbind(Y_RNAr,Y_H3K4r,Y_ATACr)



  # fit NMF model
  fit_NMF = NNLM::nnmf(X,k=10,method='scd',loss='mkl',max.iter = 50)

  init = list(L_init = fit_NMF$W,F_init = t(fit_NMF$H))
  X = Matrix(X,sparse = T)
  # run stm + nugget
  fit_stm = stm(X,K=10,nugget=TRUE,printevery = 10,tol=1e-4,init = init)
  fit_hals = NMF_HALS(as.matrix(X),K=10,smooth_method = 'runmed',printevery = 10)
  saveRDS(list(gene=gene,
               merge_len=merge_len,
               X=X,
               assays = c('RNA,H3K4me1','ATAC'),
               fit_NMF = fit_NMF,
               fit_stm = fit_stm,
               fit_hals=fit_hals),file=paste('output/luis/',gene,'_K10_merge10base.rds',sep=''))

}




