rm(list=ls()) 
options(stringsAsFactors = F) 

setwd('D:\\项目\\子宫内膜癌\\04.Suv')
dir.create("results")
gene=read.delim("Glu.txt",header = F)
gene=as.matrix(gene)
tcga_dat=read.delim("../01.Data pre/TCGA/results/tcga_dat.txt",sep='\t',header = T,check.names = F)
gene_set=list(gene)
names(gene_set)='Gc'
library(genefilter)
library(GSVA)
library(Biobase)
library(edgeR)
# gsva方法
gsva_matrix<- gsva(as.matrix(tcga_dat), gene_set,
                   method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)


gsva_matrix=as.data.frame(gsva_matrix)
## 肿瘤和正???
data=gsva_matrix
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group_list=ifelse(group=="0",'tumor','normal')
group_list=factor(group_list,levels = c('normal','tumor'))
rt=as.data.frame(t(gsva_matrix))
rt$group=group_list


library(ggpubr)

library(tidyverse)
colnames(rt)

p1=ggboxplot(rt,x = "group",
          y = "Gc",
          color = "black",
          fill = "group",
          xlab = "group",
          ylab = "Gc", palette=c('#00AFBB','#FC4E07')
)+stat_compare_means()
ggsave(filename = "results/Boxplot.pdf",p1,he=5,wi=5)
dev.off()
tcga_cli_T=read.delim("../01.Data pre/TCGA/results/tcga_cli_T.txt",sep='\t',header = T,check.names = F)
gene_set = read.table(file='marker.txt', header = T, sep = '\t',stringsAsFactors = F)
list <- list()
for(i in 1:length(unique(gene_set$cluster))){
  list[[i]] <- gene_set$gene[gene_set$cluster== (unique(gene_set$cluster)[i])]
}
names(list)<- unique(gene_set$cluster)
gsva_matrix<- gsva(as.matrix(tcga_dat), list,
                   method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)


gsva_matrix=as.data.frame(gsva_matrix)
data=gsva_matrix
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group_list=ifelse(group=="0",'tumor','normal')
group_list=factor(group_list,levels = c('normal','tumor'))

rt=as.data.frame(t(gsva_matrix))
rt$group=group_list
table(rt$group)
library(ggpubr)

library(tidyverse)
rt=tidyr::pivot_longer(rt,cols = -c('group'),names_to = "Gc",values_to = 'Abundance')

library(ggpubr)

p2=ggboxplot(
  rt,
  x = "Gc",
  y = "Abundance",
  color = "black",
  fill = "group",
  xlab = "group",
  ylab = "Abundance", palette=c('#00AFBB','#FC4E07')
) +
  stat_compare_means(
    aes(group = group),
    label = "p.signif", 
    method = "wilcox.test",
    hide.ns = T,
    size = 4
  ) +
  theme(axis.text.x = element_text(
    angle =45,
    hjust = 1,
    vjust = 1
  ))
ggsave(filename = "results/celltype_Boxplot.pdf",p2,he=5,wi=5)
dev.off()


data=gsva_matrix
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[,group==0]

data=as.data.frame(t(data))

#训练集
tcga_dat_T=read.delim("../01.Data pre/Dataset/train.data.txt",sep='\t',header = T,check.names = F)
tcga_cli_T=read.delim("../01.Data pre/Dataset/train.cli.txt",sep='\t',header = T,check.names = F)

list <- list()
for(i in 1:length(unique(gene_set$cluster))){
  list[[i]] <- gene_set$gene[gene_set$cluster== (unique(gene_set$cluster)[i])]
}
names(list)<- unique(gene_set$cluster)
gsva_matrix<- gsva(as.matrix(tcga_dat_T), list,
                   method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)


gsva_matrix=as.data.frame(gsva_matrix)
data=as.data.frame(t(gsva_matrix))

# 最优Cutoff寻找预后影响
## 读取生存数据
suv=tcga_cli_T
cli=dplyr::select(suv,'OS.time','OS')
colnames(cli)=c("futime", "fustat")
rownames(cli)=suv$sample
identical(rownames(cli),rownames(data))
range(data)
# 0.3096521 1.3096521
out=cbind(cli,data)
out=cbind(id=row.names(out),out)
write.table(out,file="results/train_expTime.txt",sep="\t",row.names=F,quote=F)


rt=read.table("results/train_expTime.txt", header=T, sep="\t", check.names=F, row.names=1)     #读取输入文件

library(survminer)
library(survival)
library(dplyr)
library(caret)
pFilter=1                                               
outTab=data.frame()
sigGenes=c("futime","fustat")

for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<pFilter){
    sigGenes=c(sigGenes,i)
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}

write.table(outTab,file="results/uniCox_train.txt",sep="\t",row.names=F,quote=F)



## K-M
immune_p=c()
immune_figure=list()
library(survival)
library(survminer)

#install.packages('cowplot')
library(cowplot)
dev.off()
dir.create('k_m')
range(rt$futime)
ggplotKM<-function(time,status,group,labs,palette,leg){
  library(ggplot2)
  library(survival)
  dat1=data.frame(time=time,status=status,group=group)
  colnames(dat1)=c('time','status','groups')
  sf<-survival::survfit(Surv(time,status) ~ groups,data=dat1)
  surv=survminer::ggsurvplot(sf, data = dat1, 
                             palette = palette, 
                             pval = TRUE,
                             surv.median.line='hv'
                             #,conf.int = T
                             ,conf.int.style ='step'
                             , pval.coord=c(0, 0.2), #Add p-value 
                             risk.table = TRUE, 
                             legend.title = leg
                             ,legend.labs =labs
                             ,conf.int=T
  )
  p1=surv$plot+theme_bw()+
    theme(axis.text.y=element_text(family="Times",face="plain"),
          axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          plot.margin=unit(c(0.2, 0.2, 0, 0.1), "inches"),
          legend.position=c(1,1),
          legend.justification=c(1,1),
          legend.background = element_rect(fill = NA, colour = NA),
          legend.title = element_text(family="Times",face="plain"),
          legend.text = element_text(family="Times",face="plain"))
  p2=surv$table+theme_bw()+
    theme(axis.text.y=element_text(family="Times",face="plain"),
          plot.margin=unit(c(0, 0.2, 0.2, 0.1), "inches"),
          plot.title=element_blank(),
          legend.position=c(1,1), 
          legend.justification=c(1,1),
          legend.title = element_text(family="Times",face="plain"),
          legend.text = element_text(family="Times",face="plain"))+
    xlab('Time(Year)')
  
  g2=ggpubr::ggarrange(p1,p2, ncol = 1, nrow = 2,heights = c(1,0.3),align = "v")
  return(g2)
}
for (i in rownames(gsva_matrix)) {
  res.cut=surv_cutpoint(rt, time="futime", event="fustat", variables=i)
  cutoff=as.numeric(res.cut$cutpoint[1])
  print(cutoff)
  Type=ifelse(rt[,i]<= cutoff, "Low", "High")
  data=rt[,c("futime","fustat")]
  data$group=Type
  data$group=factor(data$group, levels=c("Low", "High"))
  bioCol=c('#00AFBB','#FC4E07')
  p=ggplotKM(time = data$futime,status = data$fustat,group = data$group,
             labs = c('Low','High'),palette = bioCol,leg = i)
  ggsave(paste0('k_m/',i,'.pdf'),p,height = 6,width = 6)
}

#验证集
test_dat=read.delim("../01.Data pre/Dataset/test.data.txt",sep='\t',header = T,check.names = F)

list <- list()
for(i in 1:length(unique(gene_set$cluster))){
  list[[i]] <- gene_set$gene[gene_set$cluster== (unique(gene_set$cluster)[i])]
}
names(list)<- unique(gene_set$cluster)
gsva_matrix2<- gsva(as.matrix(test_dat), list,
                   method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_matrix2=as.data.frame(gsva_matrix2)


data2=as.data.frame(t(gsva_matrix2))

test_cli=read.table('../01.Data pre/Dataset/test.cli.txt',row.names = 1,header = T,check.names = F)
test_cli=dplyr::select(test_cli,'OS.time','OS')
colnames(test_cli)=c("futime", "fustat")

identical(rownames(data2),rownames(test_cli))
out2=cbind(test_cli,data2)
out2=cbind(id=row.names(out2),out2)
write.table(out2,file="results/expTime_test.txt",sep="\t",row.names=F,quote=F)

rt=read.table("results/expTime_test.txt", header=T, sep="\t", check.names=F, row.names=1)
pFilter=1                                               
outTab=data.frame()
sigGenes=c("futime","fustat")


for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<pFilter){
    sigGenes=c(sigGenes,i)
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}
write.table(outTab,file="results/uniCox_test.txt",sep="\t",row.names=F,quote=F)
rt=read.table("results/expTime_test.txt", header=T, sep="\t", check.names=F, row.names=1) 


immune_p=c()
immune_figure=list()
library(survival)
library(survminer)
library(cowplot)
dir.create('test1')
ggplotKM<-function(time,status,group,labs,palette,leg){
  library(ggplot2)
  library(survival)
  dat1=data.frame(time=time,status=status,group=group)
  colnames(dat1)=c('time','status','groups')
  sf<-survival::survfit(Surv(time,status) ~ groups,data=dat1)
  surv=survminer::ggsurvplot(sf, data = dat1, 
                             palette = palette, 
                             pval = TRUE,
                             surv.median.line='hv'
                             #,conf.int = T
                             ,conf.int.style ='step'
                             , pval.coord=c(0, 0.2), #Add p-value 
                             risk.table = TRUE, 
                             legend.title = leg
                             ,legend.labs =labs
                             ,conf.int=T
  )
  p1=surv$plot+theme_bw()+
    theme(axis.text.y=element_text(family="Times",face="plain"),
          axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          plot.margin=unit(c(0.2, 0.2, 0, 0.1), "inches"),
          legend.position=c(1,1),
          legend.justification=c(1,1),
          legend.background = element_rect(fill = NA, colour = NA),
          legend.title = element_text(family="Times",face="plain"),
          legend.text = element_text(family="Times",face="plain"))
  p2=surv$table+theme_bw()+
    theme(axis.text.y=element_text(family="Times",face="plain"),
          plot.margin=unit(c(0, 0.2, 0.2, 0.1), "inches"),
          plot.title=element_blank(),
          legend.position=c(1,1), 
          legend.justification=c(1,1),
          legend.title = element_text(family="Times",face="plain"),
          legend.text = element_text(family="Times",face="plain"))+
    xlab('Time(Year)')
  
  g2=ggpubr::ggarrange(p1,p2, ncol = 1, nrow = 2,heights = c(1,0.3),align = "v")
  return(g2)
}

for (i in rownames(gsva_matrix)) {
  res.cut=surv_cutpoint(rt, time="futime", event="fustat", variables=i)
  cutoff=as.numeric(res.cut$cutpoint[1])
  print(cutoff)
  Type=ifelse(rt[,i]<= cutoff, "Low", "High")
  data=rt[,c("futime","fustat")]
  data$group=Type
  data$group=factor(data$group, levels=c("Low", "High"))
  bioCol=c('#00AFBB','#FC4E07')
  p=ggplotKM(time = data$futime,status = data$fustat,group = data$group,
             labs = c('Low','High'),palette = bioCol,leg = i)
  ggsave(paste0('test/',i,'.pdf'),p,height = 6,width = 6)
}


#UNICOX
rt1=read.table('results/uniCox_train.txt',header = T)
rt2=read.table('results/uniCox_test.txt',header = T)
#rt3=read.table('results/uniCox_METABRIC.txt',header = T)
rt1$Dataset='Train'
rt2$Dataset='Test'
#rt3$Dataset='METABRIC'


rt=rbind(rt1,rt2)
rt$`-logP`=-log(rt$pvalue)
rt$logHR=log(rt$HR)
library(ggplot2)
rt$logHR=ifelse(rt$logHR>3,3,rt$logHR)
rt$logHR=ifelse(rt$logHR< -3,-3,rt$logHR)


ggplot(data=rt)+
  geom_point(aes(y=Dataset,x=id,fill=logHR,size= `-logP`),
             color='black',shape=21,stroke=1.5)+
  scale_fill_gradientn(colours = c('#403D76','#E3B635','#C02E20'),limits=c(-3,3))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=60,hjust=1),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank())+
  scale_size_area(breaks=c(1,2,3))
