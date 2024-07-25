Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())
setwd('D:\\项目\\子宫内膜癌/07.Single gene immune/')
dir.create('results')
hubgene=read.delim("../06.Immune/modelgene.txt",header = F)
filter_genes=hubgene$V1
train_exp <- read.delim("../01.Data pre/Dataset/train.data.txt",sep='\t',header = T,check.names = F,row.names = 1)
train_cli=read.delim("../01.Data pre/Dataset/train.cli.txt",sep='\t',header = T,check.names = F,row.names = 1)
#使用训练集进行分析
library(psych)
library(clusterProfiler)
hallmark.list <- clusterProfiler::read.gmt('msigdb_v7.5_GMTs/h.all.v7.5.symbols.gmt')
hallmark.list <- split(hallmark.list$gene, hallmark.list$term)
ssGSEA_score<-function(dat, kcdf = c("Gaussian", "Poisson")[1], GeneSet =  NULL) {
  # dat：行为基因，列为样本
  # kcdf：log2(FPKM、TPM)数据使用 Gaussian；Counts 数据使用 Poisson
  library(GSVA)
  
  if (is.null(GeneSet)) {
    gene_set <- read.table(paste0(basedir,"/SourceFiles/PMID_28052254.txt"), header = T, stringsAsFactors = F, sep = '\t')[, 1:2]
    gene_set <- split(as.matrix(gene_set)[,1], gene_set[,2])
  } else {
    gene_set <- GeneSet
  }
  
  ssgsea_res <- gsva(as.matrix(dat), 
                     gene_set,
                     method='ssgsea',
                     kcdf=kcdf,
                     abs.ranking=TRUE)
  ssgsea_res <- as.data.frame(t(ssgsea_res))
  return(ssgsea_res)
}


train_hallmark <- ssGSEA_score(dat = train_exp,
                               GeneSet = hallmark.list)
dim(train_hallmark)
colnames(train_hallmark) <- gsub('HALLMARK_', '', colnames(train_hallmark))
train_hallmark_cor <- corr.test(t(train_exp[filter_genes, ]),
                                as.matrix(train_hallmark[colnames(train_exp), ]),
                                method="pearson",
                                adjust="none", 
                                ci=F)
train_hallmark_cor$r
train_hallmark_cor_R <- train_hallmark_cor$r
train_hallmark_cor_P <- train_hallmark_cor$p
library(gtools)
train_hallmark_cor_P  <- stars.pval(train_hallmark_cor_P)
colnames(train_hallmark_cor_P) <- colnames(train_hallmark_cor_R)
# install.packages('superheat')
library(superheat)
pdf('results/pathways.cor.res.pdf', width = 10, height = 10, onefile = F)
mycolor <- ggsci::pal_igv(alpha =1)(8)
superheat(t(train_hallmark_cor_R),
          X.text=t(train_hallmark_cor_P),
          bottom.label.text.angle = 0,
          pretty.order.rows = TRUE,
          pretty.order.cols = TRUE,
          bottom.label.text.size = 3,
          left.label.text.size = 3,
          left.label.size = 0.5,
          bottom.label.size = 0.05,
          heat.pal = colorRampPalette(c(rep("navy",1), "white", rep("firebrick3",1)))(100))
dev.off()


# 与免疫评分的相关性分析
load("wb_funtions.Rdata")
train_estimate <- estimate_score(dat = train_exp, platform = "affymetrix")
head(train_estimate)
library(ggstatsplot)
for (ge in filter_genes) {
  dir.create(paste0('results/gene2estimate/', ge), 
             recursive = T, showWarnings = F)
  for (ic in colnames(train_estimate)[1:3]){
    tmp.data <- data.frame(X = as.numeric(train_exp[ge, ]),
                           Y = as.numeric(train_estimate[colnames(train_exp), ic]),
                           stringsAsFactors = F)
    
    tmp.plot <- ggscatterstats(data = tmp.data, 
                               x = X,
                               y = Y, 
                               margins = "both", 
                               xfill = mycolor[1], 
                               yfill = mycolor[2], 
                               xlab = paste0('Gene expression level(', ge, ')'),
                               ylab = ic,
                               marginal = F,
                               ggtheme = theme_bw(),
                               title = NULL)
    ggsave(plot = tmp.plot,
           filename = paste0('results/gene2estimate/', 
                             ge, '/', ge, '-', ic, '.pdf'),
           width = 6, height = 5)
    ggsave(plot = tmp.plot,
           filename = paste0('results/gene2estimate/', 
                             ge, '/', ge, '-', ic, '.jpeg'),
           width = 6, height = 5)
  }
}

library(ggplot2)

MCPcounter_score<-function(dat, featuresType = "HUGO_symbols") {
  # dat：行为基因，列为样本
  res <- MCPcounter::MCPcounter.estimate(dat,
                                         featuresType = featuresType,
                                         probesets=read.table("immu_mcp_probes.txt",
                                                              sep="\t",
                                                              stringsAsFactors=FALSE,
                                                              colClasses="character"),
                                         genes=read.table("immu_mcp_genes.txt",
                                                          sep="\t",
                                                          stringsAsFactors=FALSE,
                                                          header=TRUE,
                                                          colClasses="character",
                                                          check.names=FALSE))
  res <- as.data.frame(t(res))
  return(res)
}
train_mcpcounter <- MCPcounter_score(dat = train_exp)
head(train_mcpcounter)
for (ge in filter_genes) {
  dir.create(paste0('results/gene2mcpcounter/', ge), 
             recursive = T, showWarnings = F)
  for (ic in colnames(train_mcpcounter)){
    tmp.data <- data.frame(X = as.numeric(train_exp[ge, ]),
                           Y = as.numeric(train_mcpcounter[colnames(train_exp), ic]),
                           stringsAsFactors = F)
    
    tmp.plot <- ggscatterstats(data = tmp.data, 
                               x = X,
                               y = Y, 
                               margins = "both", 
                               xfill = mycolor[1], 
                               yfill = mycolor[2], 
                               xlab = paste0('Gene expression level(', ge, ')'),
                               ylab = ic,
                               marginal = F,
                               ggtheme = theme_bw(),
                               title = NULL)
    ggsave(plot = tmp.plot,
           filename = paste0('results/gene2mcpcounter/', 
                             ge, '/', ge, '-', ic, '.pdf'),
           width = 6, height = 5)
    ggsave(plot = tmp.plot,
           filename = paste0('results/gene2mcpcounter/', 
                             ge, '/', ge, '-', ic, '.jpeg'),
           width = 6, height = 5)
  }
}


ssGSEA_score<-function(dat, kcdf = c("Gaussian", "Poisson")[1], GeneSet =  NULL) {
  # dat：行为基因，列为样本
  # kcdf：log2(FPKM、TPM)数据使用 Gaussian；Counts 数据使用 Poisson
  library(GSVA)
  
  if (is.null(GeneSet)) {
    gene_set <- read.table("PMID_28052254.txt.txt", header = T, stringsAsFactors = F, sep = '\t')[, 1:2]
    gene_set <- split(as.matrix(gene_set)[,1], gene_set[,2])
  } else {
    gene_set <- GeneSet
  }
  
  ssgsea_res <- gsva(as.matrix(dat), 
                     gene_set,
                     method='ssgsea',
                     kcdf=kcdf,
                     abs.ranking=TRUE)
  ssgsea_res <- as.data.frame(t(ssgsea_res))
  return(ssgsea_res)
}

train_ssgsea_28 <- ssGSEA_score(dat = train_exp)
head(train_ssgsea_28)
for (ge in filter_genes) {
  dir.create(paste0('results/gene2PMID_28052254/', ge), 
             recursive = T, showWarnings = F)
  for (ic in colnames(train_ssgsea_28)){
    tmp.data <- data.frame(X = as.numeric(train_exp[ge, ]),
                           Y = as.numeric(train_ssgsea_28[colnames(train_exp), ic]),
                           stringsAsFactors = F)
    
    tmp.plot <- ggscatterstats(data = tmp.data, 
                               x = X,
                               y = Y, 
                               margins = "both", 
                               xfill = mycolor[1], 
                               yfill = mycolor[2], 
                               xlab = paste0('Gene expression level(', ge, ')'),
                               ylab = ic,
                               marginal = F,
                               ggtheme = theme_bw(),
                               title = NULL)
    ggsave(plot = tmp.plot,
           filename = paste0('results/gene2PMID_28052254/', 
                             ge, '/', ge, '-', ic, '.pdf'),
           width = 6, height = 5)
    ggsave(plot = tmp.plot,
           filename = paste0('results/gene2PMID_28052254/', 
                             ge, '/', ge, '-', ic, '.jpeg'),
           width = 6, height = 5)
  }
}


train_immune_score <- cbind(train_ssgsea_28,
                            # GSE75214_cibersort,
                            train_mcpcounter,
                            train_estimate[, 1:3])
dim(train_immune_score)
colnames(train_immune_score)
dim(train_mcpcounter)
train_immune_score_group <- data.frame(Features = colnames(train_immune_score),
                                       Groups = c(rep('ssGSEA', 28),
                                                  # rep('CIBERSORT', 22),
                                                  rep('MCP-Counter', 10),
                                                  rep('ESTIMATE', 3)),
                                       stringsAsFactors = F)

dim(train_immune_score)
train_immune_score_cor <- corr.test(t(train_exp[filter_genes, ]),
                                    as.matrix(train_immune_score[colnames(train_exp), ]),
                                    method="pearson",
                                    adjust="none", 
                                    ci=F)

train_immune_score_cor_R <- round(train_immune_score_cor$r, 2)

train_immune_score_cor_P <- train_immune_score_cor$p
train_immune_score_cor_P <- ifelse(train_immune_score_cor_P < 0.0001, '****',
                                   ifelse(train_immune_score_cor_P < 0.001, '***',
                                          ifelse(train_immune_score_cor_P < 0.01, "**",
                                                 ifelse(train_immune_score_cor_P < 0.05, '*', 'ns'))))


train_immune_score_cor_P1 <- paste0(train_immune_score_cor_R, '\n', train_immune_score_cor_P)

train_immune_score_cor_P1 <- matrix(train_immune_score_cor_P1, nrow = 4)
train_immune_score_cor_P1=train_immune_score_cor_P1[,]
colnames(train_immune_score_cor_P1) <- colnames(train_immune_score_cor_P)
rownames(train_immune_score_cor_P1) <- rownames(train_immune_score_cor_P)


library(ComplexHeatmap)
mycolors=c("#5050FFFF","#CE3D32FF","#749B58FF","#F0E685FF","#466983FF","#BA6338FF","#5DB1DDFF",
           "#802268FF","#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F",
           "#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928","#1B9E77","#D95F02","#7570B3",
           "#E7298A","#66A61E","#E6AB02","#A6761D","#666666","#E41A1C","#377EB8","#4DAF4A",
           "#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999")

pdf('results/gene_immune_corr.pdf', width = 15, height = 15)
Heatmap(t(train_immune_score_cor_R),
        name = "Corr",
        col = circlize::colorRamp2(c(-1, 0, 1),
                                   c(mycolors[1], 'white', mycolors[2])),
        border = T,
        show_column_names = F,
        show_column_dend = F,
        show_row_dend = F,
        cluster_columns=T,
        cluster_rows=T,
        column_names_side = 'bottom',
        column_names_rot = 60,
        row_names_side = 'right',
        column_split = factor(rownames(train_immune_score_cor_R)),
        row_split = factor(train_immune_score_group$Groups),
        row_title_rot = 0,
        top_annotation = HeatmapAnnotation(Gene = rownames(train_immune_score_cor_R)
                                           , col=list(Gene=c('MTHFD2'=mycolors[1],
                                                             'MRPS18C' = mycolors[2],
                                                             'YBX1'=mycolors[3],
                                                             'ANXA2'=mycolors[4]
                                                             ))
                                           , annotation_width = unit(c(1,2), 'cm')
                                           , simple_anno_size = unit(1, "cm")
                                           , annotation_height = unit(0.2, "cm")
                                           , gap = unit(1, 'mm')),
        left_annotation = rowAnnotation(Groups = train_immune_score_group$Groups
                                        , col=list(Groups=c('ESTIMATE'=mycolors[6],
                                                            'MCP-Counter'=mycolors[7],
                                                            'ssGSEA'=mycolors[8]))
                                        , annotation_width = unit(c(1,2), 'cm')
                                        , annotation_height = unit(0.2, "cm")
                                        , simple_anno_size = unit(3, "cm")
                                        , gap = unit(1, 'mm')),
        # rect_gp = gpar(col = NA, lwd = 1),
        row_names_gp = gpar(fontsize = 12),
        rect_gp = gpar(col = "white", lwd = 2),
        cell_fun = function(j, i, x, y, w, h, col) {
          grid.text(t(train_immune_score_cor_P1)[i, j], x, y)
        }
)
dev.off()




#筛选重要分子
library(org.Hs.eg.db)
library(GOSemSim)
library(reshape2)
library(ggplot2)
entriz <- mapIds(org.Hs.eg.db, keys = filter_genes, keytype = "SYMBOL", column="ENTREZID")
entriz 
rt=data.frame(ENTREZID=entriz,SYMBOL=filter_genes)
rt$ENTREZID <- as.character(rt$ENTREZID)
head(rt)
bp <- godata('org.Hs.eg.db', ont="BP", computeIC = FALSE)
cc <- godata('org.Hs.eg.db', ont="CC", computeIC = FALSE)
mf <- godata('org.Hs.eg.db', ont="MF", computeIC = FALSE)
simbp <- mgeneSim(rt$ENTREZID,
                  semData = bp,
                  measure = "Wang",
                  drop = NULL,
                  combine = "BMA")
simcc <- mgeneSim(rt$ENTREZID,
                  semData = cc,
                  measure = "Wang",
                  drop = NULL,
                  combine = "BMA")
simmf <- mgeneSim(rt$ENTREZID,
                  semData = mf,
                  measure = "Wang",
                  drop = NULL,
                  combine = "BMA")
fsim <- (simmf * simcc * simbp)^(1/3)
colnames(fsim) = rt$SYMBOL
rownames(fsim) = rt$SYMBOL
for (i in 1:ncol(fsim)){
  fsim[i,i] <- NA
}
dat <- melt(fsim)
dat <- dat[!is.na(dat$value),]
dat <- dat[,c(1,3)]
head(dat)
dat.mean <- aggregate(value~Var1, dat, mean)

m <- dat.mean$value
names(m) <- dat.mean$Var1
#按平均值给基因名排序
dat$Var1 <- factor(dat$Var1,
                   levels=names(sort(m)))
str(dat)
ggplot(dat,
       aes(x = Var1, y = value, fill = factor(Var1))) + 
  scale_fill_brewer(palette="Set2") +   #配色选择
  geom_boxplot() +
  coord_flip() +   #坐标轴互换
  xlab("") + ylab("") +
  theme_bw()

