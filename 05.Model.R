rm(list=ls()) 
options(stringsAsFactors = F) 
setwd('D:\\项目\\子宫内膜癌\\05.Model')
marker = read.table(file='../04.Suv/marker.txt', header = T, sep = '\t',stringsAsFactors = F)
marker=marker[1:180,]
diff=read.delim("../03.Gene scelect/results/deg.txt",header = T)
#diff=diff[which(diff$adj.P.Val<0.01),]
comgene=intersect(marker$gene,rownames(diff))                
#abs(diff$logFC)>0.5&diff$
#得到43个差异基因，但基因数量仍然太多，我们需要进行cox差异分析得到更佳的预后基
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
write.table(tcga.sig.fit,file = "results/modelgene.txt",quote = F,sep='\t',row.names = T,col.names=T)


modelgene=rownames(tcga.sig.fit)
#训练集
traindata=read.delim("../01.Data pre/Dataset/train.data.txt",check.names = F)
traindata=traindata[modelgene,]
traincli=read.delim("../01.Data pre/Dataset/train.cli.txt",header = T)

write.table(traindata,file = "results/traindata.txt",quote = F,sep='\t',col.names = T)
write.table(traincli,file = "results/traincli.txt",quote = F,sep='\t',col.names = T,row.names = F)

#验证集
testdata=read.delim("../01.Data pre/Dataset/test.data.txt",check.names = F)
testdata=testdata[modelgene,]
testcli=read.delim("../01.Data pre/Dataset/test.cli.txt",header = T)

write.table(testdata,file = "results/testdata.txt",quote = F,sep='\t',col.names = T)
write.table(testcli,file = "results/testcli.txt",quote = F,sep='\t',col.names = T,row.names = F)



#计算出model risk
model_risk=read.delim("./Model/riskscore_mat.txt",header = T,row.names = 1,check.names = F)
model_risk=model_risk[,"RSF",drop=F]
#先是训练集
train_risk=model_risk[traincli$sample,,drop=F]
train_risk$group=ifelse(train_risk$RSF>median(train_risk$RSF),"High","Low")

train_risk=cbind(traincli,train_risk)


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
  p1=surv$plot+theme_minimal()+
    theme(axis.text.y=element_text(family="Times",face="plain"),
          axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          plot.margin=unit(c(0.2, 0.2, 0, 0.1), "inches"),
          legend.position=c(1,1),
          legend.justification=c(1,1),
          legend.background = element_rect(fill = NA, colour = NA),
          legend.title = element_text(family="Times",face="plain"),
          legend.text = element_text(family="Times",face="plain"))
  p2=surv$table+theme_minimal()+
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
#提取肿瘤临床信息，看一下高低干性组的预后
library(survminer) 
library(survival) 
#仍然使用res.cat计算
risk_plot<-function(time,event,riskscore,
                    group,mk,labs=c('High','Low'),
                    palette=ggsci::pal_aaas()(10)){
  dat=data.frame(time=time,status=event,riskscore=riskscore,group=group)
  #ROC曲线
  ROC_rt=timeROC::timeROC(T=dat$time, 
                          delta=dat$status,
                          marker=dat$riskscore, cause=1,
                          weighting='marginal',
                          times=mk, 
                          ROC=TRUE,iid = T)
  p.dat=data.frame()
  for(i in which(ROC_rt$times>0)){
    lbs=paste0('AUC at ',mk[i],' years: ',round(ROC_rt$AUC[i],2))
    p.dat=rbind.data.frame(p.dat,data.frame(V1=ROC_rt$FP[,i],V2=ROC_rt$TP[,i],Type=lbs))
  }
  #colnames(p.dat)=c('V1','V2','Type')
  p.dat=as.data.frame(p.dat)
  p.dat$Type=as.factor(p.dat$Type)
  roc_plot=ggplot(p.dat, aes(x=V1,y=V2, fill=Type))+
    stat_smooth(aes(colour=Type),se = FALSE, size = 1)+
    theme_bw()+xlab('False positive fraction')+
    ylab('True positive fraction') +
    theme(axis.text.y=element_text(family="Times",face="plain"),
          axis.text.x=element_text(family="Times",face="plain"),
          axis.title.x=element_text(family="Times",face="plain"),
          axis.title.y=element_text(family="Times",face="plain"),
          plot.title=element_blank(),
          plot.margin=unit(c(0.2, 0.2, 0.2, 0.2), "inches"),
          legend.position=c(1,0),
          legend.justification=c(1,0),
          legend.background = element_rect(fill = NA, colour = NA),
          legend.title = element_text(family="Times",face="plain"),
          legend.text = element_text(family="Times",face="plain"))
  #KM曲线
  library(ggplot2)
  library(survival)
  dat1=dat[,c("time","status","group")]
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
  gg<-ggpubr::ggarrange(g2,roc_plot,ncol = 2,nrow = 1,labels = c('A','B'))
  return(gg)
}
train_roc_km<-risk_plot(time=train_risk$OS.time,
                       event=train_risk$OS,
                       riskscore=train_risk$RSF,
                       group=train_risk$group,
                       mk=c(1,3,5),labs=c('High','Low'),
                       palette=c("#FC4E07","#00AFBB"))
train_roc_km
ggsave(filename = "results/train_roc_km.pdf",train_roc_km,he=4,wi=8)
dev.off()

rt=train_risk
library(survival)
library(survminer)
rt$OS.time=as.numeric(rt$OS.time)
res.cut=surv_cutpoint(rt, time="OS.time", 
                      event="OS",
                      variables='RSF')
cutoff=as.numeric(res.cut$cutpoint[1])
print(cutoff)
Type=ifelse(rt[,'RSF']<= cutoff, "Low", "High")
data=rt
data$group=Type
data$group=factor(data$group, levels=c("Low", "High"))
diff=survdiff(Surv(OS.time, OS) ~ group, data = data)
length=length(levels(factor(data[,"group"])))
pValue=1-pchisq(diff$chisq, df=length-1)
fit <- survfit(Surv(OS.time, OS) ~ group, data = data)
bioCol=c("#00AFBB","#FC4E07")
bioCol=bioCol[1:length]
ggsurvplot(fit, 
           data=data,
           conf.int=T,
           pval=pValue,
           pval.size=6,
           legend.title='Riskscore',
           legend.labs=levels(factor(data[,"group"])),
           legend = c(0.7, 0.8),
           font.legend=12,
           xlab="Time(Years)",
           palette = bioCol,
           #surv.median.line = "hv",
           risk.table=T,
           cumevents=F,
           risk.table.height=.25)






#验证集
test_risk=model_risk[testcli$sample,,drop=F]
test_risk$group=ifelse(test_risk$RSF>median(test_risk$RSF),"High","Low")

test_risk=cbind(testcli,test_risk)
test_roc_km<-risk_plot(time=test_risk$OS.time,
                        event=test_risk$OS,
                        riskscore=test_risk$RSF,
                        group=test_risk$group,
                        mk=c(1,3,5),labs=c('High','Low'),
                        palette=c("#FC4E07","#00AFBB"))
test_roc_km
ggsave(filename = "results/test_roc_km.pdf",test_roc_km,he=4,wi=8)
dev.off()
