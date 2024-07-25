rm(list=ls()) 
options(stringsAsFactors = F) 
setwd('D:\\项目\\子宫内膜癌\\01.Data pre\\TCGA')
dir.create('results')
tcga_cli<-read.delim('TCGA-UCEC.survival.tsv',sep='\t',header = T)

tcga_cli$OS.time=as.numeric(tcga_cli$OS.time)
tcga_cli$OS.time=tcga_cli$OS.time/12
tcga_cli=tcga_cli[,-3]
table(tcga_cli$OS)
tcga_cli=unique(tcga_cli)
rownames(tcga_cli)=tcga_cli$sample


#表达谱
tcga_dat<-read.delim('UCEC_TCGA.merge.mRNA.TPM.txt',sep='\t',header = T,check.names = F)
tcga_dat[1:4,1:4]
table(tcga_dat$type)
tcga_dat=tcga_dat[which(tcga_dat$type=='protein_coding'),]
tcga_dat=tcga_dat[,-c(1,3)]
tcga_dat=aggregate(.~Symbol,tcga_dat,mean)
rownames(tcga_dat)=tcga_dat$Symbol
tcga_dat$Symbol=NULL


tcga_cli$sample=substr(tcga_cli$sample,1,15)
colnames(tcga_dat)=substr(colnames(tcga_dat),1,15)


sam_T<-colnames(tcga_dat)[substr(colnames(tcga_dat),14,15)=='01']
sam_T=intersect(sam_T,tcga_cli$sample)
sam_N<-colnames(tcga_dat)[substr(colnames(tcga_dat),14,15)=='11']
tcga_dat_T=tcga_dat[,sam_T]
tcga_dat_N=tcga_dat[,sam_N]
tcga_cli=tcga_cli[sam_T,]
tcga_cli=na.omit(tcga_cli)
rownames(tcga_cli)=tcga_cli$sample
dim(tcga_dat)
range(tcga_dat)
tcga_dat=log2(tcga_dat+1)
#提取肿瘤样本
length(sam_T)
tcga_dat_T=tcga_dat_T[,tcga_cli$sample]

identical(colnames(tcga_dat_T),rownames(tcga_cli))
write.table(tcga_dat_T,file = 'results/tcga_dat_T.txt',quote = F,sep='\t',row.names = T)
write.table(tcga_cli,file = 'results/tcga_cli_T.txt',quote = F,sep='\t',row.names = F)

tcga_dat=tcga_dat[,c(sam_N,colnames(tcga_dat_T))]
group=as.data.frame(colnames(tcga_dat))
group$group=ifelse(substr(colnames(tcga_dat),14,15)=='01',"Tumor","Normal")
colnames(group)=c("Samples","group")
write.table(tcga_dat,file = 'results/tcga_dat.txt',quote = F,sep='\t',row.names = T)
write.table(group,file = 'results/group.txt',quote = F,sep='\t',row.names = F)




