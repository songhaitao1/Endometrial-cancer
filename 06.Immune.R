Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
setwd("D:/项目/子宫内膜癌/06.Immune")
dir.create("results")
library(limma) 
library(dplyr)
library(tidyverse)
#设定输入的文件及路径
expFile='D:/项目/子宫内膜癌/01.Data pre/Dataset/train.data.txt'
dir.create('results')
dir.create('results\\estimate')
#estimate评分
library(utils)
library(estimate)
tcga_exp=read.delim('D:/项目/子宫内膜癌/01.Data pre/Dataset/train.data.txt',sep='\t',header = T,check.names = F)
filterCommonGenes(input.f=expFile , output.f="results\\estimate\\estimate_input.gct", id="GeneSymbol")#过滤掉细胞中共有的基因
### 计算免疫浸润打分 ###
estimateScore("results\\estimate\\estimate_input.gct", "results\\estimate\\estimate_score.gct", platform="affymetrix")#计算肿瘤微环境的评分
ESTIMATE_score = read.table("results\\estimate\\estimate_score.gct",
                            skip = 2,
                            header = T,
                            row.names = 1,check.names = F)
ESTIMATE_score = as.data.frame(t(ESTIMATE_score[,2:ncol(ESTIMATE_score)]))
ESTIMATE_score$Samples = rownames(ESTIMATE_score)
ESTIMATE_score = ESTIMATE_score[,c(ncol(ESTIMATE_score),2:ncol(ESTIMATE_score)-1)]
row.names(ESTIMATE_score)=colnames(tcga_exp)
ESTIMATE_score=ESTIMATE_score[,-1]
write.csv(ESTIMATE_score, "results\\estimate\\ESTIMATE_score.CSV")



#计算出model risk
model_risk=read.delim("../05.Model/Model/riskscore_mat.txt",header = T,row.names = 1,check.names = F)
model_risk=model_risk[,"RSF",drop=F]
#先是训练集
train_risk=model_risk[colnames(tcga_exp),,drop=F]
train_risk$group=ifelse(train_risk$RSF>median(train_risk$RSF),"High","Low")

ESTIMATE_score=ESTIMATE_score[row.names(train_risk),]
ESTIMATE_score$cluster=train_risk$group
str(ESTIMATE_score)
library(ggpubr)
library(ggsci)
head(ESTIMATE_score)
my_comparisons=list(c("High", "Low"))
PurityScore=ESTIMATE_score %>%
  ggplot(aes(x=cluster, y=TumorPurity, fill=cluster))+xlab(NULL) + ylab("PurityScore")+
  geom_violin(trim=FALSE,color="white",width = 0.6) +
  scale_x_discrete(labels = c("High",'Low'))+
  geom_boxplot(width=0.2,position=position_dodge(0.9))+scale_fill_manual(values = c("#FC4E07","#00AFBB"))+theme_bw()+
  theme_bw()+geom_blank(aes(y = 1.3))+theme(legend.position = 'none')+theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+
  stat_compare_means(comparisons=my_comparisons,method="wilcox.test",label = "p.signif")
PurityScore
ggsave(PurityScore,filename = "./results/PurityScore.pdf",he=6,wi=5)

ImmuneScore=ESTIMATE_score %>%
  ggplot(aes(x=cluster, y=ImmuneScore, fill=cluster))+xlab(NULL) + ylab("ImmuneScore")+
  geom_violin(trim=FALSE,color="white",width = 0.6) +
  scale_x_discrete(labels = c("High",'Low'))+
  geom_boxplot(width=0.2,position=position_dodge(0.9))+scale_fill_manual(values = c("#FC4E07","#00AFBB"))+theme_bw()+
  theme_bw()+geom_blank(aes(y = 1.3))+theme(legend.position = 'none')+theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+
  stat_compare_means(comparisons=my_comparisons,method="wilcox.test",label = "p.signif")
ImmuneScore
ggsave(ImmuneScore,filename = "./results/ImmuneScore.pdf",he=5,wi=5)

ESTIMATEscore=ESTIMATE_score %>%
  ggplot(aes(x=cluster, y=ESTIMATEScore, fill=cluster))+xlab(NULL) + ylab("ESTIMATEScore")+
  geom_violin(trim=FALSE,color="white",width = 0.6) +
  scale_x_discrete(labels = c("High",'Low'))+
  geom_boxplot(width=0.2,position=position_dodge(0.9))+scale_fill_manual(values = c("#FC4E07","#00AFBB"))+theme_bw()+
  theme_bw()+geom_blank(aes(y = 1.3))+theme(legend.position = 'none')+theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+
  stat_compare_means(comparisons=my_comparisons,method="wilcox.test",label = "p.signif")
ESTIMATEscore
ggsave(ESTIMATEscore,filename = "./results/ESTIMATEscore.pdf",he=5,wi=5)

StromalScore=ESTIMATE_score %>%
  ggplot(aes(x=cluster, y=StromalScore, fill=cluster))+xlab(NULL) + ylab("StromalScore")+
  geom_violin(trim=FALSE,color="white",width = 0.6) +
  scale_x_discrete(labels = c("High",'Low'))+
  geom_boxplot(width=0.2,position=position_dodge(0.9))+scale_fill_manual(values = c("#FC4E07","#00AFBB"))+theme_bw()+
  theme_bw()+geom_blank(aes(y = 1.3))+theme(legend.position = 'none')+theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+
  stat_compare_means(comparisons=my_comparisons,method="wilcox.test",label = "p.signif")
StromalScore
ggsave(StromalScore,filename = "./results/StromalScore.pdf",he=5,wi=5)




##下面用IOBR包绘制多种免疫浸润差异复杂热图
gc()
library(IOBR)
#第一步仍然是读入表达矩阵: 行名为人的基因symbol，列为样本
eset=read.table(expFile, header=T, sep="\t", check.names=F)

eset<-as.matrix(eset)
eset[1:5, 1:5]
####################################
array<-FALSE  #是否是芯片的数据，如果是RNAseq, 此处为FALSE
ProjectID<-"TCGA-UCEC"
tumor_type<-"UCEC" #命名癌种用于TIMER的解析
#各种微环境工具的解析,这里是所有相关免疫微环境计算方法，选择需要的进行计算
cibersort<-deconvo_tme(eset = eset,method = "cibersort",arrays = array,perm = 1000 )
epic     <-deconvo_tme(eset = eset,method = "epic",arrays = array)
mcp      <-deconvo_tme(eset = eset,method = "mcpcounter")
timer    <-deconvo_tme(eset = eset,method = "timer",group_list = rep(tumor_type,dim(eset)[2]))
tme_combine<-cibersort %>% 
  inner_join(.,mcp,by       = "ID") %>% 
  inner_join(.,epic,by      = "ID") %>% 
  inner_join(.,timer,by     = "ID")
# %>%inner_join(.,xcell,by     = "ID") %>%  inner_join(.,estimate,by  = "ID") %>%   inner_join(.,quantiseq,by = "ID") %>% inner_join(.,ips,by= "ID")

write.table(tme_combine,file="results/tme_combine.txt",sep="\t",quote=F,col.names=T)       

#构建细胞方法数据框
tme_combine=read.table("results/tme_combine.txt",row.names=1,check.names = F,header = T)
cibersort_method=colnames(cibersort[2:ncol(cibersort)])
cibersort_method
cibersort_method=as.data.frame(cibersort_method)
cibersort_method$method="cibersort"
colnames(cibersort_method)=c("Name", "methods")
head(cibersort_method)

epic_method=colnames(epic[2:ncol(epic)])
epic_method
epic_method=as.data.frame(epic_method)
epic_method$method="EPIC"
colnames(epic_method)=c("Name", "methods")
head(epic_method)

mcp_method=colnames(mcp[2:ncol(mcp)])
mcp_method
mcp_method=as.data.frame(mcp_method)
mcp_method$method="MCPcounter"
colnames(mcp_method)=c("Name", "methods")
head(mcp_method)

timer_method=colnames(timer[2:ncol(timer)])
timer_method
timer_method=as.data.frame(timer_method)
timer_method$method="TIMER"
colnames(timer_method)=c("Name", "methods")
head(timer_method)

library(RColorBrewer)
library(circlize)
library(gplots)
library(viridis)
library(oompaBase)
#将几种方式合并
tcga.immu.methods=rbind(cibersort_method,epic_method,mcp_method,timer_method)
methods.col <- brewer.pal(n = length(unique(tcga.immu.methods$methods)),name = "Paired")
# 创建注释
# 列注释，位于热图顶端
#读取分组信息

annCol <- data.frame(Cluster = train_risk$group,
                     row.names = rownames(train_risk),
                     stringsAsFactors = F)
# 行注释，位于热图左侧
annRow <- data.frame(Methods = factor(tcga.immu.methods$methods,
                                      levels = unique(tcga.immu.methods$methods)),
                     row.names = tcga.immu.methods$Name,
                     stringsAsFactors = F)

standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}
# 数据标准化
tme_combine=as.data.frame(tme_combine)
rownames(tme_combine)=tme_combine$ID
tme_combine=tme_combine[,-1]
indata <- t(tme_combine)
indata <- indata[,colSums(indata) > 0] # 确保没有富集全为0的细胞
plotdata <- standarize.fun(indata,halfwidth = 2)
lab_rownames=rownames(plotdata)
num=as.numeric(table(annCol$Cluster))
diff_cell<-function(dat,group){
  dat=data.frame(cluster=group,t(dat))
  gr=as.character(unique(group))
  dat1=t(dat[dat$cluster==gr[1],-1])
  dat2=t(dat[dat$cluster==gr[2],-1])
  pathway=unique(c(rownames(dat1),rownames(dat2)))
  p_vale=data.frame()
  for (i in pathway){
    dd1=t.test(as.numeric(dat1[i,]),as.numeric(dat2[i,]))$p.value
    p_vale=rbind(p_vale,data.frame(pathway=i,p.value=dd1))
  }
  return(p_vale)
}
tme_combine=tme_combine[rownames(annCol),]
pathwy_p<-diff_cell(dat=t(tme_combine),group=annCol$Cluster)
head(pathwy_p)
pathwy_p_fit=pathwy_p[pathwy_p$p.value<0.05,]
pathwy_p_fit$lab=ifelse(pathwy_p_fit$p.value<0.001,'***',ifelse(pathwy_p_fit$p.value<0.01,'**',ifelse(pathwy_p_fit$p.value<0.05,'*','neg')))
pathwy_p_fit$value=ifelse(pathwy_p_fit$p.value<0.001,'<0.001',round(pathwy_p_fit$p.value,3))
plotdata1=as.data.frame(plotdata)

plotdata1$pathway=rownames(plotdata1)
plotdata_res=merge(data.frame(pathway=pathwy_p_fit$pathway,
                              lab=paste0(pathwy_p_fit$pathway,' ',
                                         pathwy_p_fit$lab,' ',
                                         pathwy_p_fit$value)),
                   plotdata1,by='pathway')
plotdata_res=plotdata_res[,-1]
rownames(plotdata_res)=plotdata_res$lab
plotdata_res$lab=gsub('_TIMER','',plotdata_res$lab)
plotdata_res$lab=gsub('_CIBERSORT','',plotdata_res$lab)
plotdata_res$lab=gsub('_EPIC','',plotdata_res$lab)
plotdata_res$lab=gsub('_MCPcounter','',plotdata_res$lab)

annRow=merge(annRow,pathwy_p_fit,by.x = 0,by.y ="pathway" )
head(annRow)
colnames(annRow)[1]=c("pathway")
annRow_res=merge(data.frame(pathway=annRow$pathway,lab=paste0(annRow$pathway,' ',annRow$lab,' ',annRow$value)),annRow,by='pathway')
annRow_res=annRow_res[,c(2,3)]
rownames(annRow_res)=annRow_res$lab.x
annRow_res=annRow_res[,-1,drop=F]

library(pheatmap)
anno_col=data.frame(Cluster=train_risk$group)
rownames(anno_col)=row.names(train_risk)
anno_col=anno_col[order(anno_col$Cluster),,drop=F]
colnames(annRow_res)[1]='Methods'
annRow_res=annRow_res[order(annRow_res$Methods),,drop=F]
dat=plotdata_res[rownames(annRow_res),]
table(annRow_res$Methods)
table(anno_col$Cluster)

library(ggsci)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
pdf('results/immune heatmap by pheatmap.pdf',height = 9,width = 12)
pheatmap(mat=as.matrix(dat[rownames(annRow_res),rownames(anno_col)]),
         scale = 'row',show_colnames = F,
         show_rownames = T,cluster_cols = F,annotation_row = annRow_res,
         annotation_col = anno_col,labels_row = dat$lab,cluster_rows = F,
         gaps_col = c(187),gaps_row = c(6,13,19),
         annotation_colors = list(Cluster=c('High'="#FC4E07",
                                            'Low'="#00AFBB"),
                                  Methods=c('cibersort'=pal_nejm()(8)[1],
                                            'EPIC'=pal_nejm()(8)[2],
                                            'MCPcounter'=pal_nejm()(8)[3],
                                            'TIMER'=pal_nejm()(8)[4])),
         color=c(colorRampPalette(colors=c("blue","white"))(length(bk)/2),
                 colorRampPalette(color=c("white","red"))(length(bk)/2)),
         legend_breaks=seq(-2,2,1),breaks=bk,annotation_names_row = F)
dev.off()



#1、IMV210数据集####
imv210.exp<-read.delim('IMvigor210/IMvigor210_Counts2TPM.txt',sep='\t',header = T)
IMvigor210_anno <- read.delim('IMvigor210/IMvigor210_entrez2gene.txt', header = T)
imv210.exp=merge(IMvigor210_anno[,c("entrez_id","symbol")],
                 imv210.exp,by.x='entrez_id',by.y='genes')
imv210.exp[1:4,1:4]
imv210.exp=imv210.exp[,-1]
imv210.exp=aggregate(.~symbol,imv210.exp,mean)
imv210.exp[1:4,1:4]
rownames(imv210.exp)=imv210.exp$symbol
imv210.exp=imv210.exp[-1,-1]
range(imv210.exp)
imv210.exp <- log2(imv210.exp + 1)
IMvigor210_cli <- read.delim('IMvigor210/IMvigor210_cli.txt', header = T)
IMvigor210_cli <- IMvigor210_cli[, c("os", "censOS",
                                     "Best.Confirmed.Overall.Response",
                                     "IC.Level", "TC.Level", "Immune.phenotype",
                                     "FMOne.mutation.burden.per.MB",
                                     "Neoantigen.burden.per.MB","TCGA.Subtype")]
colnames(IMvigor210_cli) <- c('OS.time', 'OS', 'Response', 
                              "IC.Level", "TC.Level", "Immune.phenotype",
                              'TMB', 'NEO','Stage')
table(IMvigor210_cli$Response)
IMvigor210_cli=IMvigor210_cli[which(IMvigor210_cli$Response != 'NE'),]
IMvigor210_cli$Response[IMvigor210_cli$Response=='CR'|IMvigor210_cli$Response=='PR']<-'CR/PR'
IMvigor210_cli$Response[IMvigor210_cli$Response=='PD'|IMvigor210_cli$Response=='SD']<-'PD/SD'
imv210.exp=imv210.exp[,rownames(IMvigor210_cli)]
module.gene=read.delim("modelgene.txt",header = F)
module.gene=module.gene$V1
get_riskscore<-function(dat,os,os.time,step=T,direction=c("both", "backward", "forward")[1]){
  crbind2DataFrame=function(dat){
    print(class(dat))
    if(class(dat)=='table'){
      if(!is.na(ncol(dat))){
        dat=apply(dat,2,function(x){
          return(x)
        })
      }
    }
    if(class(dat)!='data.frame'){
      dat1=as.data.frame(as.matrix(dat))
    }else{
      dat1=dat
    }
    #print(head(dat1))
    for(i in 1:ncol(dat1)){
      dat1[,i]=as.character(dat1[,i])
      dt=dat1[which(gsub(' ','',dat1[,i])!=''&!is.na(dat1[,i])),i]
      dt=dt[which(dt!='Inf'&dt!='NaN'&dt!='NA')]
      if(sum(is.na(as.numeric(dt)))==0){
        #print(dat1[,i])
        dat1[,i]=as.numeric(dat1[,i])
      }
    }
    return(dat1)  
  }
  
  tcga_dat1 <- cbind(time=os.time,
                     status=os,
                     dat)
  tcga_dat1=crbind2DataFrame(tcga_dat1)
  colnames(tcga_dat1)=gsub('-','__',colnames(tcga_dat1))
  gene111=gsub('-','__',colnames(dat))
  fmla <- as.formula(paste0("Surv(time, status) ~"
                            ,paste0(gene111,collapse = '+')))
  cox <- coxph(fmla, data =as.data.frame(tcga_dat1))
  if(step==T){
    cox1 <- step(cox,direction =direction)
  }else{
    cox1=cox
  }
  lan <- coef(cox1)
  #round(lan, 3)
  genes <- names(cox1$coefficients)
  mult_results=paste0(round(lan, 3), '*', names(lan),collapse = '+')
  risk.tcga <- as.numeric(lan%*%as.matrix(t(tcga_dat1[,genes])))
  
  data_gene_score_final<-tcga_dat1
  data_gene_score_final$Samples<-rownames(data_gene_score_final)
  data_gene_score_final$riskscore=risk.tcga
  data_gene_score_final$riskscorez=mosaic::zscore(risk.tcga)
  optimalCutoff <- survminer::surv_cutpoint(data.frame(time=data_gene_score_final$time/365,
                                                       event=data_gene_score_final$status,
                                                       risk=data_gene_score_final$riskscore), 
                                            time = "time", event = "event",variables = c("risk"))
  optimalCutoff=optimalCutoff$cutpoint$cutpoint[1]
  #print(optimalCutoff)
  #optimalCutoff=median(data_gene_score_final$riskscore)
  #optimalCutoff=0
  data_gene_score_final$Risk=ifelse(data_gene_score_final$riskscore>optimalCutoff,'High','Low')
  table(data_gene_score_final$Risk)
  data_gene_score_final$cutoff=optimalCutoff
  return(list(result=data_gene_score_final,module.gene=cox1$coefficients,model=mult_results))
}
imv.risk<-get_riskscore(dat = as.data.frame(t(imv210.exp[intersect(rownames(imv210.exp),module.gene),])),
                        os = IMvigor210_cli$OS,
                        os.time = IMvigor210_cli$OS.time,
                        step = F,direction = 'both')
imv.risk$result$Risk=ifelse(imv.risk$result$riskscorez>0,'High','Low')
library(ggplot2)
library(survival)
ggplotKM<-function(time,status,group,labs,palette){
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
                             legend.title = 'Group'
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
          legend.text = element_text(family="Times",face="plain"))
  
  g2=ggpubr::ggarrange(p1,p2, ncol = 1, nrow = 2,heights = c(1,0.3),align = "v")
  return(g2)
}
fig1a<-ggplotKM(time = imv.risk$result$time/365,
                status =imv.risk$result$status,
                group = imv.risk$result$Risk,
                labs = c('High','Low'),
                palette = c("#FC4E07","#00AFBB"))
fig1a
dev.off()
ggsave(fig1a,filename = "./results/imkm.pdf",he=5,wi=5)


imv.risk.cli<-merge(imv.risk$result[,c("Samples","riskscorez","Risk")],
                    data.frame(Samples=rownames(IMvigor210_cli),IMvigor210_cli),
                    by='Samples')
write.table(x = imv.risk.cli,'results/imv.risk.cli.txt',quote = F,row.names = F,sep='\t')
Response=imv.risk.cli[,c("Response","riskscorez")]
response.risk=prop.table(table(imv.risk.cli$Response,imv.risk.cli$Risk),margin=2)
table(Response$Response)
my_comparisons=list(c("CR/PR",'PD/SD'))
Response=Response %>%
  ggplot(aes(x=Response, y=riskscorez, fill=Response))+xlab(NULL) + ylab("Response")+
  geom_violin(trim=FALSE,color="white",width = 0.6) +
  scale_x_discrete(labels = c("CR/PR",'PD/SD'))+
  geom_boxplot(width=0.2,position=position_dodge(0.9))+scale_fill_manual(values = c("#FC4E07","#00AFBB"))+theme_bw()+
  theme_bw()+geom_blank(aes(y = 1.3))+theme(legend.position = 'none')+theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+
  stat_compare_means(comparisons=my_comparisons,method="wilcox.test",label = "p.signif")
Response
ggsave(Response,filename = "./results/Response.pdf",he=5,wi=5)
dev.off()


response.risk=reshape2::melt(response.risk)
colnames(response.risk)<-c("type","Risk","Percentage")

response.risk$Percentage<-round(response.risk$Percentage,digits=2)
write.table(response.risk,'results/response.risk.txt',quote = F,row.names = F,sep='\t')
pdf('./results/Response比例.pdf')
fig1c=ggplot(response.risk,aes(x=Risk,y=Percentage,fill=type))+
  geom_bar(position = "fill",stat="identity")+
  theme_bw()+
  geom_text(aes(label = Percentage),position=position_stack(vjust =0.5),size = 5)+
  ggsci::scale_fill_nejm()
fig1c
dev.off()
imv.risk.cli2=imv.risk.cli[which(imv.risk.cli$Stage=='III'|imv.risk.cli$Stage=='IV'),]
pdf('./results/IMV III-IV.pdf',wi=5,he=5)
fig1e<-ggplotKM(time = imv.risk.cli2$OS.time/30,
                status =imv.risk.cli2$OS,
                group = imv.risk.cli2$Risk,
                labs = c('High','Low'),
                palette = c("#FC4E07","#00AFBB"))#定义颜色
fig1e
dev.off()





#TIP
table(train_risk$group)
cluster=train_risk[,2,drop=F]
colnames(cluster)="cluster"
TIP=read.delim("./TIP/ssGSEA.normalized.score.txt",sep='\t',header = T,check.names = F,row.names = 1)
#原文中主要聚焦于T细胞，我们以Step4.T cell.recruiting为后续分析
TIP=TIP[c(1:4,21:23),]
rownames(TIP)=c("Step1","Step2","Step3","Step4","Step5","Step6","Step7")
TIP=as.data.frame(t(TIP))
TIP=TIP[row.names(cluster),]
TIP=cbind(TIP,cluster)

#绘制雷达图
library(tidyverse) 
library(ggradar)  
library(palmerpenguins)
library(scales)


radar <- TIP %>%  #管道符
  drop_na() %>%  #去掉有NA值的行
  group_by(cluster) %>%  #根据不同species分组，用于后面summarise()计算
  summarise(
    Step1 = mean(Step1),
    Step2 = mean(Step2),
    Step3 = mean(Step3),
    Step4 = mean(Step4),
    Step5 = mean(Step5),
    Step6 = mean(Step6),
    Step7 = mean(Step7)
  ) %>%
  ungroup() %>%  #取消分组
  mutate_at(vars(-cluster),rescale)  #scales包的rescale()函数标准化数据到[0,1]范围
radar

T.cell=ggradar(
  radar, 
  values.radar = c("0", "0.5", "1.0"),
  grid.min = 0, grid.mid = 0.5, grid.max = 1.0,
  # Polygons
  group.line.width = 1, 
  group.point.size = 3,
  group.colours = c("#FC4E07", "#00AFBB"),
  # Background and grid lines
  background.circle.colour = "white",
  gridline.mid.colour = "grey"
)+theme_bw() +
  theme(plot.title = element_text(size = 15),
        axis.text = element_blank(),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        legend.position = "right")
T.cell
dev.off()
ggsave(T.cell,filename = "./results/T.cell.pdf",wi=5,he=4)




#这里是为了submap
traindata=read.delim("../01.Data pre/Dataset/train.data.txt",check.names = F)

cluster=cluster[colnames(traindata),,drop=F]
traindata=traindata[,rownames(cluster)]
identical(colnames(traindata),rownames(cluster))
write.table(traindata,'submap/testdata.txt',quote = F,row.names = T,sep='\t')
write.table(cluster,'submap/test.group.txt',quote = F,row.names = T,sep='\t')



#读取数据画热图
Bonferroni=t(read.delim("submap output/submap_nominal.p.matrix.Fisher.txt",sep='\t',header = T))
Normal=t(read.delim("submap output/submap_nominal.p.matrix.ES.B.on.A.txt",sep='\t',header = T))
submap=rbind(Bonferroni,Normal)
rownames(submap)=c("High","Low"," High"," Low")

colnames(submap)=c("PD1-CR","PD1-PD","PD1-PR","PD1-SD")

submap

library(ComplexHeatmap)
pdf("results/figure9c.pdf",width = 8,height = 6)
pheatmap(submap, 
         border_color = "white",
         number_format = "%.3f",
         cellwidth = 30, cellheight = 30,
         cluster_rows = F,cluster_cols = F,
         display_numbers = T,
         number_color = "black",
         fontsize_number = 9,
         name = "Statitic",
         annotation_row = data.frame(pvalue=c("Bonferroni corrected",
                                              "Bonferroni corrected",
                                              "Normal p value","Normal p value"),
                                     row.names = rownames(submap)),
         annotation_colors = list(pvalue=c("Bonferroni corrected"="#E6E6FA","Normal p value"="grey50")))

dev.off()

