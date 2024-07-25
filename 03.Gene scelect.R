rm(list=ls()) 
options(stringsAsFactors = F) 
dir.create("results")
setwd('D:\\项目\\子宫内膜癌\\03.Gene scelect')
library(dplyr)
library(tidyverse)
gene=read.delim("GT.txt",header = T)
gene=gene$GT
tcga_dat=read.delim("D:\\项目\\子宫内膜癌\\01.Data pre\\TCGA\\results\\tcga_dat.txt",sep='\t',header = T,check.names = F)
group=read.delim("D:\\项目\\子宫内膜癌\\01.Data pre\\TCGA\\results\\group.txt",sep='\t',header = T,check.names = F)
table(group$group)
group_list=c(rep('Normal',35),rep('Tumor',536))
## 强制限定顺序
group_list <- factor(group_list,levels = c("Normal","Tumor"),ordered = F)
exprSet <- tcga_dat
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
#差异分析
dat <- exprSet
design=model.matrix(~factor(group_list))
fit=lmFit(dat,design)
fit=eBayes(fit)
options(digits = 4)
topTable(fit,coef=2,adjust='BH')
bp=function(g){
  library(ggpubr)
  df=data.frame(gene=g,stage=group_list)
  p <- ggboxplot(df, x = "stage", y = "gene",
                 color = "stage", palette = "jco",
                 add = "jitter")
  #  Add p-value
  p + stat_compare_means()
}
deg=topTable(fit,coef=2,adjust='BH',number = Inf)
head(deg) 

#筛选logFC大于1，adj.P.value小于0.05
deg=deg[which(deg$adj.P.Val<0.05),]
write.table(deg,file = "results/deg.txt",quote = F,sep='\t',row.names = T,col.names=T)

comgene=intersect(gene,rownames(deg))


#得到111个差异基因，但基因数量仍然太多，我们需要进行cox差异分析得到更佳的预后基
tcga_dat_T=read.delim("D:\\项目\\子宫内膜癌\\01.Data pre\\TCGA\\results\\tcga_dat_T.txt",sep='\t',header = T,check.names = F)
tcga_cli_T=read.delim('D:\\项目\\子宫内膜癌\\01.Data pre\\TCGA\\results\\tcga_cli_T.txt',sep='\t',header = T)
tcga_dat_m=as.data.frame(t(tcga_dat_T[comgene,]))
tcga_dat_m=cbind.data.frame(OS.time=tcga_cli_T$OS.time,
                            OS=tcga_cli_T$OS,
                            tcga_dat_m[tcga_cli_T$sample,])

library(survival)
cox_batch<-function(dat,time,event){
  coxRun<-function(dat){
    library(survival)
    colnames(dat)=c('time','status','AS')  
    dat=dat[which(!is.na(dat[,1])&!is.na(dat[,3])&!is.na(dat[,2])),]
    #print(nrow(dat))
    if(nrow(dat)<10){
      print(paste0('Sample Num is small:',nrow(dat)))
      return(c(NA,NA,NA,NA))
    }
    #if(quantile(dat[,3])['25%']==quantile(dat[,3])['50%']) return(c(NA,NA,NA,NA))
    fmla <- as.formula("Surv(time, status) ~AS")
    if(table(dat[,2])[1]>1&table(dat[,2])[2]>1){
      cox <- survival::coxph(fmla, data = dat)
      re=c(summary(cox)[[7]][5],summary(cox)[[7]][2],summary(cox)[[8]][3],summary(cox)[[8]][4])
      return(re)
    }else{
      return(c(NA,NA,NA,NA))
    }
  }
  t.inds=which(!is.na(time)&!is.na(event))
  dat1=dat[,t.inds]
  os=time[t.inds]
  ev=event[t.inds]
  
  ct=sum(ev%in%c(0,1))
  if(ct!=length(ev)){
    print('event must be 0(alive) or 1(dead)')
    return(NULL)
  }
  
  res=t(apply(dat1, 1, function(x){
    ep=as.numeric(as.character(x))
    ind2=which(!is.na(ep))
    print(length(ind2))
    if(length(ind2)>1){
      os1=os[ind2]
      ev1=ev[ind2]
      ep1=ep[ind2]
      return(coxRun(data.frame(os1,ev1,ep1)))
    }else{
      return(c(NA,NA,NA,NA))
    }
  }))
  colnames(res)=c('p.value','HR','Low 95%CI','High 95%CI')
  row.names(res)=row.names(dat1)
  return(as.data.frame(res))
}
tcga.cox=cox_batch(dat = t(tcga_dat_m[,-c(1,2)]),
                   time = tcga_dat_m$OS.time,
                   event = tcga_dat_m$OS)
head(tcga.cox)
tcga.cox=as.data.frame(tcga.cox)
fitp=0.05
tcga.sig.fit=tcga.cox[which(tcga.cox$p.value<fitp),]

write.table(tcga.sig.fit,file = "results/tcga.sig.fit.txt",quote = F,sep='\t',row.names = T,col.names=T)
