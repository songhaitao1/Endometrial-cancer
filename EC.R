setwd("~/EC/")

options(stringsAsFactors = F)
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(ggsci)
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(stringr)
library(viridis)
library(scCustomize)
dir.create("results")
dir_name=list.dirs('GSE173682/',full.names = F,recursive = F)
dir_name
datalist=list()
#读取10x的数据创建CreateSeuratObject对象
for (i in 1:length(dir_name)){
  dir.10x = paste0("GSE173682/",dir_name[i])
  my.data <- Read10X(data.dir = dir.10x) 
  #细胞增加标签
  colnames(my.data)=paste0(dir_name[i],colnames(my.data))
  datalist[[i]]=CreateSeuratObject(counts = my.data, project = dir_name[i], min.cells = 3, min.features = 250)
  datalist[[i]]$Samples=dir_name[i]
  datalist[[i]]$type=substr(dir_name[i],1,1)
}
names(datalist)=dir_name
#批量计算线粒体和rRNA的含量
for (i in 1:length(datalist)){
  sce <- datalist[[i]]
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")# 计算线粒体占比
  sce[["percent.Ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")# 计算rRNA占比
  datalist[[i]] <- sce
  rm(sce)
}
#合并所有的数据
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
#细胞数的统计
raw_cell=sce@meta.data
raw_count <- table(raw_cell$Samples)
raw_count
sum(raw_count)#35782
pearplot_befor<-VlnPlot(sce,group.by ='Samples', 
                        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"),
                        pt.size = 0, 
                        ncol = 4)
pearplot_befor
ggsave('results/pearplot_befor.pdf',pearplot_befor,height = 5,width = 15)
ggsave('results/pearplot_befor.jpg',pearplot_befor,height = 5,width = 15,dpi = 300)
#样本的颜色
sample_color<-pal_nejm(alpha = 0.5)(8)[1:5]
sample_color
Feature_ber1<-FeatureScatter(sce,feature1 = 'nFeature_RNA',
                             feature2 = 'nCount_RNA',
                             group.by = 'Samples',
                             cols = sample_color)
Feature_ber2<-FeatureScatter(sce,feature1 = 'percent.mt',
                             feature2 = 'nCount_RNA',
                             group.by = 'Samples',
                             cols = sample_color)
Feature_ber3<-FeatureScatter(sce,feature1 = 'percent.mt',
                             feature2 = 'nFeature_RNA',
                             group.by = 'Samples',
                             cols = sample_color)
Feature_ber1=Feature_ber1+theme(legend.position = 'none')
Feature_ber2=Feature_ber2+theme(legend.position = 'none')

Feature_ber<-ggarrange(Feature_ber1,Feature_ber2,Feature_ber3,ncol = 3,nrow = 1,widths = c(1,1,1.2))
Feature_ber
ggsave('results/Feature_cor.pdf',Feature_ber,height = 5,width = 17)
ggsave('results/Feature_cor.jpg',Feature_ber,height = 5,width = 17,dpi = 300)
#过滤
datalist <- lapply(X = datalist, FUN = function(x) {
  x<-subset(x,subset =  
              nFeature_RNA < 5000 & 
              percent.mt < 15)
})
#合并数据
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
clean_cell=sce@meta.data

clean_count <- table(clean_cell$Samples)
clean_count
sum(clean_count)#30233
pearplot_after <- VlnPlot(sce,group.by ='Samples', 
                          features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Ribo"),
                          pt.size = 0, 
                          ncol = 4)
pearplot_after
ggsave('results/pearplot_after.pdf',Feature_ber,height = 5,width = 15)
ggsave('results/pearplot_after.jpg',Feature_ber,height = 5,width = 15,dpi = 300)
#保存数据
save(datalist,file = 'datalist.RData')


#过滤指标3:过滤特定基因
# Filter MALAT1 管家基因
sce <- sce[!grepl("MALAT1", rownames(sce),ignore.case = T), ]
# Filter Mitocondrial 线粒体基因
sce <- sce[!grepl("^MT-", rownames(sce),ignore.case = T), ]
# 当然，还可以过滤更多

dim(sce) 
#细胞周期评分

sce = NormalizeData(sce)
s.genes=Seurat::cc.genes.updated.2019$s.genes
g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
sce=CellCycleScoring(object = sce, 
                              s.features = s.genes, 
                              g2m.features = g2m.genes, 
                              set.ident = TRUE)
p4=VlnPlot(sce, features = c("S.Score", "G2M.Score"), group.by = "orig.ident", 
           ncol = 2, pt.size = 0)
p4
sce@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
  theme_minimal()
dim(sce)

#降维聚类跑pca
sce <- NormalizeData(sce, 
                     normalization.method = "LogNormalize",
                     scale.factor = 1e4) 
sce <- FindVariableFeatures(sce,selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(sce)
sce <- ScaleData(sce, features = scale.genes)
sce <- RunPCA(sce, features = VariableFeatures(object = sce))
DimPlot(sce, reduction = "pca", group.by = "orig.ident")
library(harmony)
seuratObj <- RunHarmony(sce, group.by.vars = "orig.ident")
names(seuratObj@reductions)
seuratObj <- RunUMAP(seuratObj,  dims = 1:50, 
                     reduction = "harmony")
DimPlot(seuratObj,reduction = "umap",label=T ) 

sce=seuratObj
sce <- FindNeighbors(sce, reduction = "harmony",
                     dims = 1:50) 
#设置不同的分辨率，观察分群效果(选择哪一个？)
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)) {
  sce=FindClusters(sce, #graph.name = "CCA_snn", 
                       resolution = res, algorithm = 1)
}
colnames(sce@meta.data)
apply(sce@meta.data[,grep("RNA_snn",colnames(sce@meta.data))],2,table)
p1_dim=plot_grid(ncol = 3, DimPlot(sce, reduction = "umap", group.by = "RNA_snn_res.0.01") + 
                   ggtitle("louvain_0.01"), DimPlot(sce, reduction = "umap", group.by = "RNA_snn_res.0.1") + 
                   ggtitle("louvain_0.1"), DimPlot(sce, reduction = "umap", group.by = "RNA_snn_res.0.2") + 
                   ggtitle("louvain_0.2"))
p1_dim
p1_dim=plot_grid(ncol = 3, DimPlot(sce, reduction = "umap", group.by = "RNA_snn_res.0.8") + 
                   ggtitle("louvain_0.8"), DimPlot(sce, reduction = "umap", group.by = "RNA_snn_res.1") + 
                   ggtitle("louvain_1"), DimPlot(sce, reduction = "umap", group.by = "RNA_snn_res.0.3") + 
                   ggtitle("louvain_0.3"))
p1_dim
library(clustree)
p2_tree=clustree(sce@meta.data, prefix = "RNA_snn_res.")
ggsave(plot=p2_tree, filename="results/Tree_diff_resolution.pdf")
#接下来分析，按照分辨率为0.8进行 
sel.clust = "RNA_snn_res.0.8"
sce <- SetIdent(sce, value = sel.clust)
table(sce@active.ident) 
saveRDS(sce, "sce.harmony.rds")
sce=readRDS("sce.harmony.rds")

DimPlot(sce,reduction = "umap",label = T)
library(ggplot2) 
genes_to_check = c('PTPRC', 'CD2','CD3D','CD3E','CD3G',"CD8A","CD4","GNLY",#T C
                   'IGKC','IGHG1','IGHA1',#B
                   'LYZ','CD68','CD14','C1QA','C1QB', #MAC
                   'VCAN','FCN1','S100A8',  # MONO
                   'MCAM','BGN','NOTCH3','GUCY1A2',#Smooth muscle cell
                   "CAPS",#Ciliated cells
                   'S100A9', 'S100A8', 'MMP19',# monocyte
                   'CLDN5','PECAM1', 'VWF','TM4SF1',  ## endo 
                   'LUM','DCN','COL6A3','FGF7','MME', 'ACTA2',#fibroblast
                   'EPCAM' , 'TACSTD2','MUC16','ALDH1A1' # epi
)
#T cell('CD2','CD3D','CD3E','CD3G'),CD4 T ('CD3D','CD4'),CD8 T('CD3D','CD8A','CD8B','GZMA'),Macrophage('CD163','CD68','CD14'),Monocyte('VCAN'),B cell('CD19','CD79A','MS4A1'),Plasma cell('CD79A','JSRP1'),Epithelial cell(EPCAM),Fibroblast('ACTA2','PDGFRB','NOTCH3'),Endothelial cell(PECAM1)

library(stringr)  
genes_to_check=str_to_upper(genes_to_check)
genes_to_check
p <- DotPlot(sce, features = unique(genes_to_check),
             assay='RNA',group.by = "RNA_snn_res.0.8")  + coord_flip()

p 

table(sce$RNA_snn_res.0.8)
sce$seurat_clusters=sce$RNA_snn_res.0.8
table(sce$seurat_clusters)


celltype=data.frame(ClusterID=0:20 ,
                    celltype= 0:20) 

celltype[celltype$ClusterID %in% c(4,8,14,20),2]='T cells'  
celltype[celltype$ClusterID %in% c(18),2]='B cells'  
celltype[celltype$ClusterID %in% c(6,15),2]='Macrophages' 
celltype[celltype$ClusterID %in% c(0,1,3,10,11,17),2]='Fibroblasts'   
celltype[celltype$ClusterID %in% c(7,12,16),2]='Endothelial' 
celltype[celltype$ClusterID %in% c(2,9,19),2]='epithelial' 
celltype[celltype$ClusterID %in% c(13),2]='Ciliated cell' 
celltype[celltype$ClusterID %in% c(5),2]='smooth muscles' 
head(celltype)
celltype
table(celltype$celltype)
sce@meta.data$celltype = "NA"
table(sce$celltype)

for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$RNA_snn_res.0.8 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce@meta.data$celltype)

Logfc = 0.5
Minpct = 0.35
DefaultAssay(sce) <- "RNA"
Idents(sce)<-'celltype'
sce.markers <- FindAllMarkers(object = sce,logfc.threshold = Logfc, min.pct = Minpct,only.pos = T)
sce.markers["pct.diff"]=sce.markers$pct.1-sce.markers$pct.2
sce.markers <- sce.markers[sce.markers$p_val_adj<0.05,]
length(unique(sce.markers$gene))
head(sce.markers)
write.table(sce.markers,'results/scRNA_marker_gene.txt',quote = F,row.names = F,sep='\t')

library(ClusterGVis)
library(org.Hs.eg.db)
pbmc.markers <- sce.markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 20, wt = avg_log2FC)
head(pbmc.markers)

st.data <- prepareDataFromscRNA(object = sce,
                                diffData = pbmc.markers,
                                showAverage = TRUE)
# enrich for clusters
library(clusterProfiler)
enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 5)
markGenes = unique(pbmc.markers$gene)[sample(1:length(unique(pbmc.markers$gene)),40,
                                             replace = F)]

# line plot
visCluster(object = st.data,
           plot.type = "line")
pdf('results/sc1.pdf',height = 10,width = 6,onefile = F)
visCluster(object = st.data,
           plot.type = "heatmap",
           column_names_rot = 45,
           markGenes = markGenes,
           cluster.order = c(1:8))
dev.off()


pdf('results/sc2.pdf',height = 10,width = 14,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich,
           line.side = "left",
           cluster.order = c(1:8),
           go.col = rep(jjAnno::useMyCol("stallion",n = 8),each = 5),
           add.bar = T)
dev.off()


#sce可视化
library(paletteer)
library(ggplot2)
pdf('results/umap.pdf',height = 4,width = 6.5,onefile = F)
pal <- c("darkgoldenrod2","deepskyblue2","#FA9FB5", "#c2e699", "#B40F20","seagreen","brown2","steelblue",
                 "#5B1A18","#9C964A","#FD6467","#6bAEd6","#fff500",
                 "#DD3497", "#7A0177", "#006837","#bcbddc", "#4a1486", "#969696","#636363","black","red")

DimPlot(sce, label = T,  cols= pal,group.by = "celltype" ,pt.size = 1, repel = T)+labs(x = "UMAP1", y = "UMAP2") +
  theme(panel.border = element_rect(fill=NA,color="black", element_rect(1), linetype="solid")) 
dev.off()
saveRDS(sce,"sce.all.RDS")
sce=readRDS("sce.all.RDS")


scRNA_T_epi=subset(sce,celltype %in% c('epithelial',"B cells"))
table(scRNA_T_epi$orig.ident)
dim(scRNA_T_epi)
epi=scRNA_T_epi
library(ggplot2)
library(Seurat)
library(dplyr)
library(reticulate)
library(sctransform)
library(cowplot)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readxl)
library(stringr)
library(cowplot)
library(scales)
library(tibble)
library(gplots)
library(RColorBrewer)

epi <- RunPCA(epi)
ElbowPlot(epi,  ndims = 50)


epi <- RunUMAP(epi, dims = 1:20)
epi <- FindNeighbors(epi, dims = 1:20)
for (i in c(0.2, 0.3, 0.4, 0.5, 1, 2)) {
  epi <- FindClusters(epi, resolution = i)
  print(DimPlot(epi, reduction = "umap", label = T) + labs(title = paste0("resolution: ", i)))
}
table(epi@meta.data$RNA_snn_res.2)
Idents(epi) <- epi@meta.data$RNA_snn_res.2
epi$cell_id=rownames(epi@meta.data)

meta_tbl <- as_tibble(FetchData(epi, c("cell_id","celltype","orig.ident")))
colnames(epi@meta.data)
pids <- sort(unique(epi$orig.ident))
pids_with_tumor <- FetchData(epi, c("celltype", "orig.ident")) %>% as_tibble() %>% 
  distinct(celltype, orig.ident) %>% filter(celltype == "B cells") %>% pull(orig.ident) %>% sort

rcm <- epi@assays$RNA@counts
dim(rcm)
af <- FetchData(epi, c("celltype", "orig.ident", "cell_id")) 

setwd("./infercnv")
## 2.准备细胞注释文件(annotation file ,简称af)
## 结果文件中，我们的malignant指的是肿瘤组织取的样
## normal指的是正常组织取的样
## 理论上，肿瘤组织取得样品有癌上皮和正常上皮两种可能
write.table(af, file = "./annotation.txt", sep = "\t", col.names = F, row.names = T)
### 3.准备基因坐标文件gene order file ,gof
### 下载网站:https://data.broadinstitute.org/Trinity/CTAT/cnv/
library(readr)
gof <- read_tsv("./gencode_v21_gen_pos.complete.txt", col_names = c("gene", "chr", "start", "end"))
### 处理下，以前有|
gof <- separate(gof,col='gene',into ='gene',sep = '\\|')
gof <- gof[!duplicated(gof$gene),]
common_genes <- intersect(gof$gene, rownames(rcm))
### 把原来的换掉
write.table(gof,file ='./gencode_v21_gen_pos.complete_modi.txt',col.names = F,sep = '\t',row.names = F,quote = F)
## 4. 指定参考分组，选择normal为参考组！tumor组相对于normal组比较（实际上是相对值）
baseline_cell_types <- "B cells"

library(dendextend)
library(infercnv)
# 构建cnv对象
cnv_obj <- infercnv::CreateInfercnvObject(
  raw_counts_matrix = rcm,
  annotations_file = "./annotation.txt",
  gene_order_file = "./gencode_v21_gen_pos.complete_modi.txt",
  ref_group_names = baseline_cell_types
) 

cnv_obj <- infercnv::run(
  cnv_obj, cutoff=0.1, out_dir="./infercnv/", 
  cluster_by_groups=T, denoise=T, HMM=T, 
  num_threads=10
)
write_rds(cnv_obj, "./infercnv/infercnv_results.rds")

## 可以加载运行结果，节约时间
cnv_obj <- readRDS("./infercnv/infercnv_results.rds")
expr <- cnv_obj@expr.data
af$cell_id=rownames(af)
af=af%>%select(cell_id,celltype)
colnames(af)=c('cell_id','group')
#因为在inferCNV过程中gof中的基因已经
gene <- rownames(expr)
gof=as.data.frame(gof)
rownames(gof)=gof$gene
# 取交集

sub_gof <-  gof[intersect(gene,gof$gene),]
expr=expr[intersect(gene,gof$gene),]
### 针对CNV推断的结果聚类
set.seed(123456)
### epi有20个cluster,此处我们也设20个好了
kmeans.result <- kmeans(t(expr), 20)
km <- data.frame(km_cluster=kmeans.result$cluster)
km$cell_id=rownames(km)

#合并km和分组信息
km=km%>%inner_join(af,by="cell_id") 
# 按照聚类结果排序
km_ordered=km[order(km$km_cluster),]
rownames(km_ordered)=km_ordered$cell_id
km_ordered$cell_id=NULL
km_ordered$km_cluster=factor(km_ordered$km_cluster) #将km_group转换为因子
head(km_ordered)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

## 准备画复杂热图
top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:22,labels_gp = gpar(cex = 1.5)))
color_v=c(RColorBrewer::brewer.pal(8, "Dark2")[1:8], RColorBrewer::brewer.pal(8, "RdGy")[1:8],
          RColorBrewer::brewer.pal(8, "Blues")[1:4])
names(color_v)=as.character(1:20)
left_anno <- rowAnnotation(df = km_ordered,col=list(group=c(
  "B cells"="#0072b5","Malignant_cells"="#bc3c27",
  'epithelial'='darkgreen'),km_cluster=color_v))


tiff("complexheatmap.tiff",width = 1200,height = 1000,units = "px")
Heatmap(t(expr)[rownames(km_ordered),], 
             col = colorRamp2(c(0.4,1,1.6), c("#377EB8","#F0F0F0","#E41A1C")), #如果是10x的数据，这里的刻度会有所变化
             cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
             column_split = factor(sub_gof$chr, paste("chr",1:22,sep = "")), #这一步可以控制染色体顺序，即使你的基因排序文件顺序是错的
             column_gap = unit(2, "mm"),
             
             heatmap_legend_param = list(title = "Modified expression",direction = "vertical",title_position = "leftcenter-rot",at=c(0.4,1,1.6),legend_height = unit(3, "cm")),
             
             top_annotation = top_anno,left_annotation = left_anno, #添加注释
             row_title = NULL,column_title = NULL)

dev.off()

#每一类对应的CB保存在kmeans_df_s数据框中
write.table(kmeans_df_s, file = "kmeans_df_s.txt", quote = FALSE, sep = '\t', row.names = T, col.names = T)
   
#我们认为5,17,20为正常细胞
table(km_ordered$km_cluster,km_ordered$group)
km_ordered$celltype_cnv=ifelse(km_ordered$km_cluster %in% c(5,17,20) & km_ordered$group==c("B cells"),"B cells",ifelse(km_ordered$km_cluster %in% c(5,17,20) & km_ordered$group==c("epithelial"),"Epithelial","Malignant"))
table(km_ordered$km_cluster,km_ordered$celltype_cnv)

metadata=epi@meta.data
cnv_type=km_ordered[,c(2,3)]

table(metadata$celltype)
table(cnv_type$celltype_cnv)

metadata=subset(metadata,celltype != "B cells")
table(metadata$celltype)
cnv_type=subset(cnv_type,celltype_cnv !="B cells")
table(cnv_type$celltype_cnv)

cnv_type=cnv_type[rownames(metadata),]
identical(rownames(cnv_type),rownames(metadata))
cnv_type=cnv_type[,-1,drop=F]
metadata=cbind(metadata,cnv_type)
table(metadata$celltype_cnv)
table(sce$celltype)
metadata2=sce@meta.data
metadata2=subset(metadata2,celltype != "epithelial")
table(metadata2$celltype)

colnames(metadata)
colnames(metadata2)
metadata=metadata[,-c(18:21)]
metadata2$celltype_cnv=metadata2$celltype

metadata_all=rbind(metadata2,metadata)
table(metadata_all$celltype_cnv)
sam=intersect(rownames(sce@meta.data),rownames(metadata_all))
length(sam)
length(rownames(sce@meta.data))

metadata_all=metadata_all[rownames(sce@meta.data),]
identical(rownames(metadata_all),rownames(sce@meta.data))

sce@meta.data=metadata_all
table(sce$celltype_cnv)
sce$celltype=sce$celltype_cnv


#打分
table(sce$celltype_cnv)
gene=read.table('Glu.txt')
gene=gene$V1
pal <- c("darkgoldenrod2","deepskyblue2","#FA9FB5", "#c2e699", "#B40F20","seagreen","brown2","steelblue",
                         "#5B1A18","#9C964A","#FD6467","#6bAEd6","#fff500",
                         "#DD3497", "#7A0177", "#006837","#bcbddc", "#4a1486", "#969696","#636363","black","red")
DoHeatmap(subset(sce,downsample=100
                 ,),features = gene,group.by = 'celltype_cnv',assay='RNA',slot = 'data',
          group.colors =pal,lines.width = 10)+
          scale_fill_gradientn(colors=c('white','firebrick3'),na.value = 'white')
metadata=sce@meta.data
dev.off()

geneset=list(gene)
names(geneset)='GTs'
gc()
library(SeuratData)
library(UCell)
library(irGSEA)
library(AUCell)
cells_rankings <- AUCell_buildRankings(sce@assays$RNA@data,  nCores=10, plotStats=TRUE) 
gc()

cells_AUC <- AUCell_calcAUC(geneset, cells_rankings,nCores =10, aucMaxRank=nrow(cells_rankings)*0.1)

aucs <- as.numeric(getAUC(cells_AUC)['GTs', ])

## 2. Ucells and singscore
sce <- irGSEA.score(object = sce, assay = "RNA",
                      slot = "data", seeds = 123, ncores = 10,msigdb = F,
                      custom = T, geneset = geneset,
                      method = c('UCell','singscore'),
                      kcdf = 'Gaussian')
uc=as.data.frame(sce@assays$UCell@counts)
uc=as.data.frame(t(uc))
singscore=as.data.frame(sce@assays$singscore@counts)
singscore=as.data.frame(t(singscore))

## 3. GSVA
exp<- as.matrix(sce@assays$RNA@data)

library(GSVA)
#开始GSVA分析
matrix = gsva(exp, 
              geneset, 
              kcdf="Gaussian",
              method="ssgsea", 
              abs.ranking=T )

ssgsea=as.data.frame(t(matrix))


# 5. addmodulescore自带函数
sce=AddModuleScore(sce,features = geneset,name = 'Add')



# 组合
score=data.frame(AUCell=aucs,UCell=uc$GTs,singscore=singscore$GTs,ssgsea=ssgsea$GTs,Add=sce$Add1)



# scale 标准化
score<- scale(score)

# 0-1标准化

normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}

score=apply(score, 2, normalize)
score=as.data.frame(score)
score$Scoring=rowSums(score)
colnames(sce@meta.data)

sce$Add1=NULL

# 将score添加入metadata，需要充分理解 
sce@meta.data=cbind(sce@meta.data,score)
save(sce,file="results/sce.all.final.RData")
#load("results/sce.all.final.RData")



library(RColorBrewer) 
library(viridis)
library(wesanderson)

n <- 30
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
pie(rep(6,n), col=sample(color, n))
col_vector
col_vector =c(wes_palette("Darjeeling1"), wes_palette("GrandBudapest1"), wes_palette("Cavalcanti1"), wes_palette("GrandBudapest2"), wes_palette("FantasticFox1"))
pal <- wes_palette("Zissou1", 12, type = "continuous")
pal2 <- wes_palette("Zissou1", 5, type = "continuous")
pal[3:12]


dev.off()
library(ggplot2)
DotPlot(sce,features = colnames(score),group.by = "celltype_cnv") + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust=1) , axis.text.y = element_text(face="bold"))+ 
  scale_colour_gradientn(colours = pal)+ theme(legend.position="right")  + labs(title = "cluster markers", y = "", x="")
dev.off()


#接下来我们聚焦于肿瘤细胞
table(sce$celltype_cnv)
sub_sce=subset(sce,celltype %in% 'Malignant')
sub_sce$id=colnames(sub_sce)
sce_Glu=sub_sce[rownames(sub_sce) %in% gene,]

df <- sce_Glu@assays$RNA@data
df=as.data.frame(df)

df <- df[rowMeans(df) !=0,  ]
df <- df[,colMeans(df) !=0 ]

## 一点没表达的细胞去除
sub_sce=subset(sub_sce,id %in% colnames(df))

library(NMF)
res <- nmf(df,10, method = "snmf/r")
sub_sce@reductions$nmf <- sub_sce@reductions$pca
sub_sce@reductions$nmf@cell.embeddings <- t(coef(res) )
sub_sce@reductions$nmf@feature.loadings <- basis(res) 
sub_sce <- RunUMAP(sub_sce,reduction = "nmf",dims = 1:10) 

sce.all <- FindNeighbors(sub_sce, reduction = "nmf",
                         dims = 1:10) 
#设置不同的分辨率，观察分群效???(选择哪一个？)
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)) {
  sce.all=FindClusters(sce.all, #graph.name = "CCA_snn", 
                       resolution = res, algorithm = 1)
}
colnames(sce.all@meta.data)
apply(sce.all@meta.data[,grep("RNA_snn",colnames(sce.all@meta.data))],2,table)
library(clustree)
clustree(sce.all@meta.data, prefix = "RNA_snn_res.")
sel.clust = "RNA_snn_res.0.5"
scRNA_T.nmf <- SetIdent(sce.all, value = sel.clust)
table(scRNA_T.nmf$RNA_snn_res.0.5)
table(Idents(scRNA_T.nmf))
## 使用nmf的分解结果降维聚???!重要
set.seed(999)

scRNA_T.nmf <- FindClusters(scRNA_T.nmf,resolution = 0.5)
scRNA_T.nmf$NMF_cluster=scRNA_T.nmf$seurat_clusters
## 结果可视???  
DimPlot(scRNA_T.nmf, label = T,cols=pal) 
                                               
df=FindAllMarkers(scRNA_T.nmf,logfc.threshold = 0.1,only.pos = T)
write.csv(df,file ='results/deg_Malignant.csv',quote=F)

head(df)
#df$avg_log2FC <- round(df$avg_log2FC, 1)
df2=df[which(df$avg_log2FC>=1),]



## 定义新的细胞亚群
NMF_celltype<- c('N-glycosylation-C1','None-C6',
                 'Elongation-C2','Other-C5','Elongation-C2','O-GalNAc-C3','O-GalNAc-C3','O-GlcNAc-C4','None-C6','None-C6','None-C6','None-C6')
Idents(scRNA_T.nmf) <- scRNA_T.nmf@meta.data$NMF_cluster
names(NMF_celltype) <- levels(scRNA_T.nmf)
scRNA_T.nmf<- RenameIdents(scRNA_T.nmf, NMF_celltype)

scRNA_T.nmf@meta.data$NMF_celltype <- Idents(scRNA_T.nmf)


Idents(scRNA_T.nmf)=scRNA_T.nmf@meta.data$NMF_celltype


DimPlot(scRNA_T.nmf,group.by = 'NMF_celltype',label = T,cols=c('#313c63','#b42e20','#ebc03e','#377b4c',
                                                                        '#7bc7cd','#5d84a4',"#8B008B"))
                                                                        ##细分
scRNA_T.nmf$celltype = factor(scRNA_T.nmf$NMF_celltype, 
                              levels=c('N-glycosylation-C1','Elongation-C2','O-GalNAc-C3','LAMTOR1-CAF-C4','O-GlcNAc-C4','Other-C5','None-C6'))


saveRDS(scRNA_T.nmf,file ='scRNA_Maligant_NMF.RDS')
scRNA_T.nmf=readRDS("scRNA_Maligant_NMF.RDS")
df=FindAllMarkers(scRNA_T.nmf,logfc.threshold = 0.25,only.pos = F)
table(df$cluster)
df$cluster = factor(df$cluster, 
                              levels=c('N-glycosylation-C1','Elongation-C2','O-GalNAc-C3','LAMTOR1-CAF-C4','O-GlcNAc-C4','Other-C5','None-C6'))
library(scRNAtoolVis)

jjVolcano(diffData = df,
          tile.col = corrplot::COL2('RdBu', 15)[4:12],
          size  = 3.5,
          fontface = 'italic',
          polar = F)
write.csv(df,file ='results/deg_Maligant_cluster.csv',quote=F)

setwd("./SCENIC/")
dim(df)
write.csv(t(as.matrix(scRNA_T.nmf@assays$RNA@counts[unique(df$gene),])),file = "sce_exp.csv")
saveRDS(sce, file = 'sce.rds')

path=getwd()

command=paste0('/home/zhangpc/anaconda3/envs/pyscenic/bin/python /home/zhangpc/pySCENIC/loom.py ',path,'/sce_exp.csv ',path)
grep_out<-system(command, intern = F)
command1=paste0('/home/zhangpc/anaconda3/envs/pyscenic/bin/pyscenic grn --num_workers 10 --sparse --method grnboost2 --output ',path,'/grn.csv ',path,'/sce.loom /home/zhangpc/pySCENIC/allTFs_hg38.txt')
grep_out<-system(command1, intern = F)
command2=paste0('/home/zhangpc/anaconda3/envs/pyscenic/bin/pyscenic ctx --num_workers 10 --output ',path,'/regulons.csv --expression_mtx_fname ',path,'/sce.loom --mode "custom_multiprocessing" --annotations_fname /home/zhangpc/pySCENIC/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl ',path,'/grn.csv ','/home/zhangpc/pySCENIC/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather')
grep_out<-system(command2, intern = F)
command3=paste0('/home/zhangpc/anaconda3/envs/pyscenic/bin/pyscenic aucell --output ',path,'/sce.sample.loom ',path,'/sce.loom ',path,'/regulons.csv --num_workers 10')
grep_out<-system(command3, intern = F)




#######
#可视化
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
#可视化相关包，多加载点没毛病
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(ggplot2)
library(pheatmap)
library(Seurat)

sce_SCENIC <- open_loom("sce.sample.loom")
regulons_incidMat <- get_regulons(sce_SCENIC, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)#将上一步矩阵文件转化为list
class(regulons)
regulonAUC <- get_regulons_AUC(sce_SCENIC, column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(sce_SCENIC)

human_data <- scRNA_T.nmf
table(human_data$NMF_celltype)
cellinfo <- human_data@meta.data[,c('NMF_celltype'),drop=F]#细胞meta信息
colnames(cellinfo)=c('celltype')
#
cellTypes <-  as.data.frame(subset(cellinfo,select = 'celltype'))
selectedResolution <- "celltype"
sub_regulonAUC <- regulonAUC
table(cellTypes)
#普通展示
next_regulonAUC <- regulonAUC[,match(colnames(human_data),colnames(regulonAUC))]
dim(next_regulonAUC)#将AUC结果于seurat对象结合

regulon_AUC <- regulonAUC@NAMES
head(next_regulonAUC[regulon_AUC,])

dd=t(next_regulonAUC@assays@data$AUC)[rownames(human_data@meta.data),regulon_AUC]

human_data@meta.data = cbind(human_data@meta.data ,dd)


TF_plot=colnames(dd)
write.table(TF_plot,'TF_list.txt',quote = F,row.names = F,sep='\t',col.names = F)

regulonActivity_byCell<-getAUC(regulonAUC)
pheatmap::pheatmap(regulonActivity_byCell[,rownames(cellinfo)], scale = 'row',annotation_col = cellinfo,
                   cluster_cols = F,show_colnames = F,
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-1, 1, length.out = 100),filename = 'regulonAUC_cell.pdf')


cellsPerGroup <- split(rownames(cellinfo), 
                       cellinfo$celltype) 

regulonActivity_byCellcluster <- sapply(cellsPerGroup,
                                        function(cells) 
                                          rowMeans(getAUC(regulonAUC)[,cells]))

regulonActivity_byCellcluster_Scaled <- t(scale(t(regulonActivity_byCellcluster), center = T, scale=T))
pheatmap::pheatmap(regulonActivity_byCellcluster_Scaled, scale = 'none',
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA,filename = 'regulonAUC_cluster.pdf')

regulons_incidMat1=t(regulons_incidMat)
regulon=data.frame()
for (i in 1:ncol(regulons_incidMat1)){
  dd=regulons_incidMat1[,i,drop=F]
  namesss=rownames(dd[which(dd[,1]>0),,drop=F])
  dd1=data.frame(TF=colnames(regulons_incidMat1)[i],gene=namesss)
  regulon=rbind.data.frame(regulon,dd1)
}

#pheatmap
write.table(regulon,'TF-gene.txt',quote = F,row.names = F,sep='\t')
# 

TF1=unique(regulon$TF)
TF1=unique(stringr::str_split_fixed(TF1,'\\(',2)[,1])
col.num=ceiling(sqrt(length(TF1)))
row.num=ceiling(length(TF1)/col.num)
#pdf('TF.umap.pdf',width = 3*col.num,height = 3*row.num)
fig1=FeaturePlot(scRNA_T.nmf,features = TF1,reduction = 'umap',pt.size = 1,ncol = col.num)+
  scale_x_continuous("")+scale_y_continuous("")+
  theme_bw()
#dev.off()
ggsave('TF.reduction.pdf',fig1,height =  3*row.num,width = 3*col.num)
ggsave('TF.reduction.jpg',fig1,height =  3*row.num,width = 3*col.num,dpi = 300)
write.table('step5.1 successful','step5.1',quote = F,row.names = F,col.names = F)



#其他可视化
#提取pyscenic第三步分析中AUC结果
regulonAUC <- get_regulons_AUC(sce_SCENIC, column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(sce_SCENIC)

##==============================加载seurat对象、RSS分析=======================================
#在可视化之前，我们再做一个分析，计算RSS值，计算regulon特异性评分
table(human_data$NMF_celltype)
cellinfo <- human_data@meta.data[,c('NMF_celltype',"nFeature_RNA","nCount_RNA")]#细胞meta信息
colnames(cellinfo)=c('NMF_celltype','nGene' ,'nUMI')
######计算细胞特异性TF
#在实际数据分析应用中，我认为比较靠谱的应用在于，细胞分了亚群，例如macrophage，有不同的特征
#我们可以查看不同亚群特异性的TF，有助于了解亚群的功能！！！！
cellTypes <-  as.data.frame(subset(cellinfo,select = 'NMF_celltype'))
selectedResolution <- "NMF_celltype"
sub_regulonAUC <- regulonAUC

rss <- calcRSS(AUC=getAUC(sub_regulonAUC),#从aucellresults获取AUC矩阵
               cellAnnotation=cellTypes[colnames(sub_regulonAUC),
                                        selectedResolution])
rss=na.omit(rss)#去除含有NA的行
#可视化细胞特异性TF
rss=rss[,c(1,4:6,2,3)]
plotRSSM_1<-function (rss, labelsToDiscard = NULL, zThreshold = 1, cluster_columns = FALSE, 
                      order_rows = TRUE, thr = 0.01, varName = "cellType", col.low = "grey90", 
                      col.mid = "darkolivegreen3", col.high = "darkgreen", revCol = FALSE, 
                      verbose = TRUE) 
{
  varSize = "RSS"
  varCol = "Z"
  if (revCol) {
    varSize = "Z"
    varCol = "RSS"
  }
  rssNorm <- scale(rss)
  rssNorm <- rssNorm[, which(!colnames(rssNorm) %in% labelsToDiscard)]
  rssNorm[rssNorm < 0] <- 0
  rssSubset <- rssNorm
  if (!is.null(zThreshold)) 
    rssSubset[rssSubset < zThreshold] <- 0
  # tmp <- .plotRSS_heatmap(rssSubset, thr = thr, cluster_columns = cluster_columns, 
  #                         order_rows = order_rows, verbose = verbose)
  # rowOrder <- rev(tmp@row_names_param$labels)
  # rm(tmp)
  rss.df <- reshape2::melt(rss)
  head(rss.df)
  colnames(rss.df) <- c("Topic", varName, "RSS")
  rssNorm.df <- reshape2::melt(rssNorm)
  colnames(rssNorm.df) <- c("Topic", varName, "Z")
  rss.df <- base::merge(rss.df, rssNorm.df)
  rss.df <- rss.df[which(!rss.df[, varName] %in% labelsToDiscard), 
  ]
  if (nrow(rss.df) < 2) 
    stop("Insufficient rows left to plot RSS.")
  #rss.df <- rss.df[which(rss.df$Topic %in% rowOrder), ]
  #rss.df[, "Topic"] <- factor(rss.df[, "Topic"], levels = rowOrder)
  p <- dotHeatmap(rss.df, var.x = varName, var.y = "Topic", 
                  var.size = varSize, min.size = 0.5, max.size = 5, var.col = varCol, 
                  col.low = col.low, col.mid = col.mid, col.high = col.high)
  invisible(list(plot = p, df = rss.df))
}

rssPlot <- 
  plotRSSM_1(rss = rss,
             zThreshold = 2,
             cluster_columns = FALSE,
             order_rows = TRUE,
             thr=0.1,
             varName = "NMF_celltype",
             col.low = '#330066',
             col.mid = '#66CC66',
             col.high = '#FFCC33')

rssPlot$plot
rssPlot$df


#提取数据，可以自己可视化dotplot，或者热图
rss_data <- rssPlot$plot$data
#devtools::install_github("XiaoLuo-boy/ggheatmap")
library(ggheatmap)
library(reshape2)
rss_data<-dcast(rss_data, 
                Topic~rss_data$NMF_celltype,
                value.var = 'Z')
rownames(rss_data) <- rss_data[,1]
rss_data <- rss_data[,-1]
colnames(rss_data)
#rss_data=rss_data[,c(1,4,5,6,2,3)]
col_ann <- data.frame(group= c(rep("N-glycosylation-C1",1),
                               rep("Elongation-C2",1),
                               rep("O-GalNAc-C3",1),
                               rep("O-GlcNAc-C4",1),
                               rep("Other-C5",1),
                               rep("None-C6",1)))#列注释
rownames(col_ann) <- colnames(rss_data)


groupcol <- c("#f57665",'#1279A2',"#CAA57D","#f1c550","#0b8457","#3B9AB2")[1:6]
names(groupcol) <- c("N-glycosylation-C1","Elongation-C2","O-GalNAc-C3","O-GlcNAc-C4","Other-C5","None-C6")
col <- list(group=groupcol)

text_columns <- sample(colnames(rss_data),0)#不显示列名
?ggheatmap
col_ann$group=factor(col_ann$group,levels = col_ann$group)
p <- ggheatmap(rss_data,color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
               cluster_rows = F,cluster_cols = F,scale = "row",
               annotation_cols = col_ann,
               annotation_color = col,levels_cols =col_ann$group,
               legendName="Relative value",
               text_show_cols = text_columns)
p
ggsave('TF.pheatmap.pdf',p,height =  8,width = 5)
dev.off()




######既然可以分析细胞中的特异性TF，那么也可以分析样本中特异性的TF
#我们可以提取某一个细胞的表达矩阵，做pyscenic，后期再进行RSS分析的时候
#可以使用样本，就可以看出哪些TF是这个样本细胞中特异性的TF，这样生物学意义也就有了
##############################################################################

##==============================TF_AUC与seurat结合===========================
#普通展示
library(Seurat)
next_regulonAUC <- regulonAUC[,match(colnames(human_data),colnames(regulonAUC))]
dim(next_regulonAUC)#将AUC结果于seurat对象结合
regulon_AUC <- regulonAUC@NAMES

human_data@meta.data = cbind(human_data@meta.data ,t(next_regulonAUC@assays@data$AUC[regulon_AUC,rownames(human_data@meta.data)]))
rownames(next_regulonAUC@assays@data$AUC[regulon_AUC,])
#自己选定感兴趣的或者比较重要的转录因子，这里我是随机的
TF_plot <- rownames(next_regulonAUC@assays@data$AUC[regulon_AUC,])

DotPlot(human_data, features = TF_plot,group.by = 'NMF_celltype')+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(hjust =1,vjust=1, angle = 45))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
####

DotPlot(human_data, features = TF_plot, group.by = 'NMF_celltype')+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(hjust =1,vjust=1, angle = 45))+
  theme(legend.direction = "horizontal", 
        legend.position = "bottom")+
  labs(x=NULL,y=NULL)
###


################################################------------------------
#------rank可视化rss-----------------------------------------------------
B_rss <- as.data.frame(rss)#rss特异性TF结果
#需要作图的细胞类型
celltype <- human_data$NMF_celltype
rssRanklist <- list()
table(celltype)
for(i in 1:length(celltype)) {
  
  data_rank_plot <- cbind(as.data.frame(rownames(B_rss)),
                          as.data.frame(B_rss[,celltype[i]]))#提取数据
  
  colnames(data_rank_plot) <- c("TF", "celltype")
  data_rank_plot=na.omit(data_rank_plot)#去除NA
  data_rank_plot <- data_rank_plot[order(data_rank_plot$celltype,decreasing=T),]#降序排列
  data_rank_plot$rank <- seq(1, nrow(data_rank_plot))#添加排序
  
  p <- ggplot(data_rank_plot, aes(x=rank, y=celltype)) + 
    geom_point(size=3, shape=16, color="#1F77B4",alpha =0.4)+
    geom_point(data = data_rank_plot[1:6,],
               size=3, color='#DC050C')+ #选择前6个标记，自行按照需求选择
    theme_bw()+
    theme(axis.title = element_text(colour = 'black', size = 12),
          axis.text = element_text(colour = 'black', size = 10),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    labs(x='Regulons Rank', y='Specificity Score',title =celltype[i])+
    ggrepel::geom_text_repel(data= data_rank_plot[1:6,],
                             aes(label=TF), color="black", size=3, fontface="italic", 
                             arrow = arrow(ends="first", length = unit(0.01, "npc")), box.padding = 0.2,
                             point.padding = 0.3, segment.color = 'black', 
                             segment.size = 0.3, force = 1, max.iter = 3e3)
  # ylim(0.3,0.6) #y轴范围
  rssRanklist[[i]] <- p
}


length(rssRanklist)
cowplot::plot_grid(rssRanklist[[1]],rssRanklist[[7]],rssRanklist[[15]],rssRanklist[[136]],rssRanklist[[2]],rssRanklist[[3]],
                   ncol=3)


#############################################################################

##==============================TF平均表达活性===========================

cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes[,selectedResolution])
#去除extend的TF
# sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 
# dim(sub_regulonAUC)

#计算平均表达
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))

#scale处理\类似于热图数据的标准化
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 


regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)

#热图
hm <- draw(ComplexHeatmap::Heatmap(regulonActivity_byGroup_Scaled, name="Regulon activity",
                                   row_names_gp=grid::gpar(fontsize=6),
                                   show_row_names = F)) 
hm
dev.off()




#=======================================================================================
#                                   pySCENIC 转录因子与靶基因
#=======================================================================================
#我们还可以可视化感兴趣的TF调控的靶基因
#TF与靶基因的关系在pyscenic分析得到的第二个文件

sce_regulons <- read.csv("regulons.csv")
sce_regulons <- sce_regulons[-2, ]
colnames(sce_regulons) <- sce_regulons[1,]
sce_regulons <- sce_regulons[-1, ]
colnames(sce_regulons) <- c("TF","ID","AUC","NES","MotifSimilarityQvalue","OrthologousIdentity",
                            "Annotation","Context","TargetGenes","RankAtMax")

#举例子我这里关注FOXP3和TGIF1这两个TF
ATF4 <- subset(sce_regulons, TF=='ATF4')
#FOXP3 <- FOXP3[which(FOXP3$AUC>0.1),]
ATF4 <- ATF4[, c("TF","TargetGenes")]
ATF4$TargetGenes <-gsub("\\[","",ATF4$TargetGenes)
ATF4$TargetGenes <-gsub("\\]","",ATF4$TargetGenes)
ATF4$TargetGenes <-gsub("\\(","",ATF4$TargetGenes)
ATF4$TargetGenes <-gsub("\\)","",ATF4$TargetGenes)
ATF4$TargetGenes <-gsub("\\'","",ATF4$TargetGenes)
library(stringr)
split_ATF4<-str_split(ATF4$TargetGenes,",")
ATF41 <- as.data.frame(split_ATF4[[1]])
# FOXP32<- as.data.frame(split_FOXP3[[2]])
# FOXP33<-as.data.frame(split_FOXP3[[3]])
# FOXP32<-as.data.frame(split_FOXP3[[4]])
# FOXP35<- as.data.frame(split_FOXP3[[5]])

names(ATF41) <- 'TF'
# names(FOXP32) <- 'TF'
# names(FOXP33) <- 'TF'
# names(FOXP34) <- 'TF'
# names(FOXP35) <- 'TF'

#FOXP3 <- rbind(FOXP31,FOXP32,FOXP33,FOXP34,FOXP35)
ATF4=ATF41

ATF4_target <- ATF4[seq(1,nrow(ATF4),2), ]
ATF4_score <- ATF4[seq(0,nrow(ATF4),2), ]

ATF4_gene <- data.frame(ATF4_target,ATF4_score)
ATF4_gene <- ATF4_gene[!duplicated(ATF4_gene$ATF4_target), ]
ATF4_gene$gene <- 'ATF4'
colnames(ATF4_gene) <- c("target","score",'tf')


FOSL1 <- subset(sce_regulons, TF=='FOSL1')
#FOXP3 <- FOXP3[which(FOXP3$AUC>0.1),]
FOSL1 <- FOSL1[, c("TF","TargetGenes")]
FOSL1$TargetGenes <-gsub("\\[","",FOSL1$TargetGenes)
FOSL1$TargetGenes <-gsub("\\]","",FOSL1$TargetGenes)
FOSL1$TargetGenes <-gsub("\\(","",FOSL1$TargetGenes)
FOSL1$TargetGenes <-gsub("\\)","",FOSL1$TargetGenes)
FOSL1$TargetGenes <-gsub("\\'","",FOSL1$TargetGenes)
library(stringr)
split_FOSL1<-str_split(FOSL1$TargetGenes,",")
FOSL11 <- as.data.frame(split_FOSL1[[1]])
# FOXP32<- as.data.frame(split_FOXP3[[2]])
# FOXP33<-as.data.frame(split_FOXP3[[3]])
# FOXP32<-as.data.frame(split_FOXP3[[4]])
# FOXP35<- as.data.frame(split_FOXP3[[5]])

names(FOSL11) <- 'TF'
# names(FOXP32) <- 'TF'
# names(FOXP33) <- 'TF'
# names(FOXP34) <- 'TF'
# names(FOXP35) <- 'TF'

#FOXP3 <- rbind(FOXP31,FOXP32,FOXP33,FOXP34,FOXP35)
FOSL1=FOSL11

FOSL1_target <- FOSL1[seq(1,nrow(FOSL1),2), ]
FOSL1_score <- ATF4[seq(0,nrow(FOSL1),2), ]

FOSL1_gene <- data.frame(FOSL1_target,FOSL1_score)
FOSL1_gene <- FOSL1_gene[!duplicated(FOSL1_gene$FOSL1_target), ]
FOSL1_gene$gene <- 'FOSL1'
colnames(FOSL1_gene) <- c("target","score",'tf')




#同理得到TGIF1及其靶基因，此处省略1万字
#two thousand years later
#TGIF1_gene
FOSL1_gene=ATF4_gene
FOSL1_gene$tf='FOSL1'
TF_target <- rbind(ATF4_gene,FOSL1_gene)
TF_target$score <- as.numeric(TF_target$score)
###接下来就是网络图了

#节点数据
paths <- c("ATF4", "FOSL1")#列重命名
nodelist <- list()
for (i in 1:length(paths)){
  node <- subset(TF_target, tf == paths[i])#提取数据
  #node=unique(node)
  nodes <- data.frame(name = unique(union(node$tf, node$target)))#整理为datafram
  nodes$value <- c(sum(node$score)/10, node$score)#加上values
  
  nodelist[[i]] <- nodes
}  #提取每个大节点数据


nodes <- rbind(nodelist[[1]],nodelist[[2]])#将三个节点文件合并
table(TF_target$tf)
table(TF_target$tf,TF_target$target)
nodes$cluster <- c(rep("ATF4",1),rep("ATF4_gene",20),
                   rep("FOSL1",1),rep("FOSL1_gene",20))#分组，为了后续颜色设置

edges <- TF_target[c("tf","target","score")]#边缘文件
edges$class <- edges$tf

library(ggraph)
library(tidygraph)
layout_cir <- tbl_graph(nodes = nodes, edges = edges)#构建ggraph作图文件
#作图
ggraph(layout_cir,layout='linear',circular = TRUE) +#选择circle
  geom_node_point(aes(size=value,colour = cluster))+#节点，大小用我们赋的值表示，颜色用分组
  geom_node_text(aes(x = 1.03 * x,
                     y = 1.03 * y,
                     label=name,
                     color=cluster,
                     angle = -((-node_angle(x, y) + 90) %% 180) + 90),
                 hjust='outward') +#文字设置。x，y是为了调整位置。angle是为了调整角度，以后其他所有网络图angle都用此公式，文字会向外发散排列
  geom_edge_arc(aes(colour=class))+#连线为曲线
  theme_void()+#theme主题
  theme(legend.position = "none")+
  scale_colour_manual(values =c('#407972',
                                         '#961E28',
                                         '#D46724',
                                         '#0f8096'))+#节点颜色
                                           scale_edge_colour_manual(values = c('#961E28',
                                                                                    '#D46724',
                                                                                        '#0f8096'))+#连线颜色
                                                                                          scale_size_continuous(range = c(2,8))+#点的大小范围设置
  coord_cartesian(xlim=c(-1.5,1.5),ylim = c(-1.5,1.5))#设置坐标位置，防止图溢出作图边界显示不全


#细胞通讯
load("results/sce.all.final.RData")
scRNA_T.nmf=readRDS("scRNA_Maligant_NMF.RDS")
table(sce$celltype_cnv)
table(scRNA_T.nmf$NMF_celltype)
sce_other=subset(sce,celltype_cnv !="Malignant")
sce_other$NMF_celltype=sce_other$celltype_cnv
table(sce_other$NMF_celltype)
colnames(sce_other@meta.data)
colnames(scRNA_T.nmf@meta.data)
scRNA_T.nmf@meta.data=scRNA_T.nmf@meta.data[,-c(25,26)]
scRNA_T.nmf@assays$UCell<- NULL
sce_other@assays$UCell<- NULL
scRNA_T.nmf@assays$singscore<- NULL
sce_other@assays$singscore<- NULL
scRNA_chat=merge(scRNA_T.nmf,sce_other)
table(scRNA_chat$NMF_celltype)
meta =scRNA_chat@meta.data # a dataframe with rownames containing cell mata data
gc()
data_input <- as.matrix(scRNA_chat@assays$RNA@data)
#data_input=data_input[,rownames(meta)]
identical(colnames(data_input),rownames(meta))

library(CellChat)
table(meta$NMF_celltype)
#O-GlcNAc-C4;LAMTOR1-CAF-C4
meta$NMF_celltype=factor(meta$NMF_celltype,levels = c('N-glycosylation-C1','Elongation-C2','O-GalNAc-C3','O-GlcNAc-C4','Other-C5','None-C6',"T cells","B cells","Ciliated cell","Endothelial","Epithelial","Fibroblasts","Myeloids","smooth muscles"))
meta=meta[order(meta$NMF_celltype),]

cellchat <- createCellChat(object = data_input[,rownames(meta)], meta = meta, group.by = "NMF_celltype")

CellChatDB <- CellChatDB.human 
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use 

dplyr::glimpse(CellChatDB$interaction)##配体-受体分析
# 提取数据库支持的数据子集
cellchat <- subsetData(cellchat)
# 识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
# 识别配体-受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
# 将配体、受体投射到PPI网络
cellchat <- projectData(cellchat, PPI.human)
unique(cellchat@idents)
cellchat@meta$NMF_celltype=as.character(cellchat@meta$NMF_celltype)
unique(cellchat@idents)
cellchat <- computeCommunProb(cellchat)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

df.net<- subsetCommunication(cellchat)

cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))


##时常deff.off!!!!
dev.off()
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= T, sources.use = c('N-glycosylation-C1','Elongation-C2','O-GalNAc-C3','O-GlcNAc-C4','Other-C5','None-C6'),
                 title.name = "Number of interactions")
dev.off()

p_bubble= netVisual_bubble(cellchat,
                           sources.use = c('N-glycosylation-C1','Elongation-C2','O-GalNAc-C3','LAMTOR1-CAF-C4','O-GlcNAc-C4','Other-C5','None-C6'),
                           remove.isolate = FALSE)+coord_flip()
p_bubble

p_bubble2= netVisual_bubble(cellchat,
                           targets.use = c('N-glycosylation-C1','Elongation-C2','O-GalNAc-C3','LAMTOR1-CAF-C4','O-GlcNAc-C4','Other-C5','None-C6'),
                           remove.isolate = FALSE)+coord_flip()
p_bubble2

dev.off()

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")


h1=netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
h2=netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
h1 + h2
dev.off()

#代谢
library(scMetabolism)
library(ggplot2)
library(rsvd)
scRNA_T.nmf$NMF_celltype <- factor(scRNA_T.nmf$NMF_celltype, levels = c('N-glycosylation-C1','Elongation-C2','O-GalNAc-C3','O-GlcNAc-C4','Other-C5','None-C6'))
scRNA_T.nmf_meta<-sc.metabolism.Seurat(obj = scRNA_T.nmf, method = "VISION", imputation = F, ncores = 2, metabolism.type = "KEGG")
input.pathway <- rownames(scRNA_T.nmf_meta@assays[["METABOLISM"]][["score"]])[1:30]
DotPlot.metabolism_1<-function (obj, pathway, phenotype, norm = "y",level.col=NULL) 
{
  input.norm = norm
  input.pathway <- pathway
  input.parameter <- phenotype
  metadata <- obj@meta.data
  metabolism.matrix <- obj@assays$METABOLISM$score
  metadata[, input.parameter] <- as.character(metadata[, input.parameter])
  metabolism.matrix_sub <- t(metabolism.matrix[input.pathway, 
  ])
  gg_table <- c()
  for (i in 1:length(input.pathway)) {
    gg_table <- rbind(gg_table, cbind(metadata[, input.parameter], 
                                      input.pathway[i], metabolism.matrix_sub[, i]))
  }
  gg_table <- data.frame(gg_table)
  gg_table_median <- c()
  input.group.x <- unique(as.character(gg_table[, 1]))
  input.group.y <- unique(as.character(gg_table[, 2]))
  for (x in 1:length(input.group.x)) {
    for (y in 1:length(input.group.y)) {
      gg_table_sub <- subset(gg_table, gg_table[, 1] == 
                               input.group.x[x] & gg_table[, 2] == input.group.y[y])
      gg_table_median <- rbind(gg_table_median, cbind(input.group.x[x], 
                                                      input.group.y[y], median(as.numeric(as.character(gg_table_sub[, 
                                                                                                                    3])))))
    }
  }
  gg_table_median <- data.frame(gg_table_median)
  gg_table_median[, 3] <- as.numeric(as.character(gg_table_median[, 
                                                                  3]))
  gg_table_median_norm <- c()
  input.group.x <- unique(as.character(gg_table[, 1]))
  input.group.y <- unique(as.character(gg_table[, 2]))
  range01 <- function(x) {
    (x - min(x))/(max(x) - min(x))
  }
  if (input.norm == "y") 
    for (y in 1:length(input.group.y)) {
      gg_table_median_sub <- subset(gg_table_median, gg_table_median[, 
                                                                     2] == input.group.y[y])
      norm_value <- range01(as.numeric(as.character(gg_table_median_sub[, 
                                                                        3])))
      gg_table_median_sub[, 3] <- norm_value
      gg_table_median_norm <- rbind(gg_table_median_norm, 
                                    gg_table_median_sub)
    }
  if (input.norm == "x") 
    for (x in 1:length(input.group.x)) {
      gg_table_median_sub <- subset(gg_table_median, gg_table_median[, 
                                                                     1] == input.group.x[x])
      norm_value <- range01(as.numeric(as.character(gg_table_median_sub[, 
                                                                        3])))
      gg_table_median_sub[, 3] <- norm_value
      gg_table_median_norm <- rbind(gg_table_median_norm, 
                                    gg_table_median_sub)
    }
  if (input.norm == "na") 
    gg_table_median_norm <- gg_table_median
  gg_table_median_norm <- data.frame(gg_table_median_norm)
  if (!is.null(level.col)){
    gg_table_median_norm[,1]=factor(gg_table_median_norm[,1],levels = level.col)
    gg_table_median_norm=gg_table_median_norm[order(gg_table_median_norm$X1),]
  }
  gg_table_median_norm[, 3] <- as.numeric(as.character(gg_table_median_norm[, 
                                                                            3]))
  library(wesanderson)
  pal <- wes_palette("Zissou1", 100, type = "continuous")
  
  ggplot(data = gg_table_median_norm, aes(x = gg_table_median_norm[, 
                                                                   1], y = gg_table_median_norm[, 2], color = gg_table_median_norm[, 
                                                                                                                                   3])) + geom_point(data = gg_table_median_norm, aes(size = gg_table_median_norm[, 
                                                                                                                                                                                                                  3])) + ylab("Metabolic Pathway") + xlab(input.parameter) + 
    theme_bw() + theme(axis.text.x = element_text(angle = 45, 
                                                  hjust = 1), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
    scale_color_gradientn(colours = pal) + labs(color = "Value", 
                                                size = "Value") + NULL
}
DotPlot.metabolism_1(obj = scRNA_T.nmf_meta,
                   pathway = input.pathway, phenotype = "NMF_celltype", norm = "y",
                   level.col= c('N-glycosylation-C1','Elongation-C2','O-GalNAc-C3','O-GlcNAc-C4','Other-C5','None-C6'))

#富集分析
library(ggforce)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
col = c("N-glycosylation-C1"="#f57665",
        "Elongation-C2"='#1279A2',
        "O-GalNAc-C3"="#CAA57D",
        "LAMTOR1-CAF-C4"="#f1c550",
        "Other-C5"="#0b8457",
        "None-C6"="#3B9AB2")

p1 = DimPlot(scRNA_T.nmf, label = T, cols = col)+NoLegend()+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())
dev.off()
#marker基因，并作富集分析，当然了，每个celltype的差异基因也是可以这样做的
sce_marker <- FindAllMarkers(scRNA_T.nmf,logfc.threshold = 0.25,
                             min.pct = 0.25,only.pos = T)

#每组celltype的基因富集分析
group <- data.frame(gene=sce_marker$gene,group=sce_marker$cluster)
Gene_ID <- bitr(sce_marker$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
data  <- merge(Gene_ID,group,by.x='SYMBOL',by.y='gene')
#GO
data_GO <- compareCluster(ENTREZID~group, 
                          data=data, 
                          fun="enrichGO", 
                          OrgDb="org.Hs.eg.db",
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.1,
                          qvalueCutoff = 0.1)

data_GO_sim <- simplify(data_GO, cutoff=0.7, by="p.adjust", select_fun=min)
data_GO_sim <- data_GO_sim@compareClusterResult
#每种celltype我们选个10条terms可视化
data_GO_sec <- data_GO_sim %>% group_by(group) %>% top_n(5, -p.adjust)
data_GO_sec$`-log10(p)` <- -log10(data_GO_sec$p.adjust)
#将需要展示的通路展示在一起
p2 = ggplot(data_GO_sec) +
  ggforce::geom_link(data = data_GO_sec, 
                     aes(x = 0,
                         y = reorder(Description, -`-log10(p)`),
                         xend =`-log10(p)`,
                         yend = Description,
                         colour=group,
                         alpha = after_stat(index),
                         size = after_stat(index)),
                     show.legend = TRUE)+
  geom_point(data = data_GO_sec,
             aes(x = `-log10(p)`, y = Description),
             size=6,
             shape = 21,
             fill = "white")+
  scale_colour_manual(values = col, 
                      guide = guide_legend(title = "Group"))+
  geom_text(mapping=aes(x=`-log10(p)`,
                        y=Description,
                        label=Count),
            size =2.5,color="black")+
  geom_vline(xintercept = -log10(0.05),
             color = "grey60",
             linetype=2, cex=0.5)+
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(colour = 'black'))

p2
dev.off()

#
#第二种富集方法
df_sig  <- sce_marker[sce_marker$p_val_adj < 0.05, ]

df_sig$cluster=factor(df_sig$cluster,levels=c('N-glycosylation-C1','Elongation-C2','O-GalNAc-C3','O-GlcNAc-C4','Other-C5','None-C6'))
table(df_sig$cluster)
library(clusterProfiler)
library(ggplot2)
group <- data.frame(gene=df_sig$gene,
                    group=df_sig$cluster)
table(df_sig$cluster)

Gene_ID <- bitr(df_sig$gene, fromType="SYMBOL", 
                toType="ENTREZID", 
                OrgDb="org.Hs.eg.db")

#构建文件并分析
data  <- merge(Gene_ID,group,by.x='SYMBOL',by.y='gene')
data_GO <- compareCluster(
  ENTREZID~group, 
  data=data, 
  fun="enrichGO", 
  OrgDb="org.Hs.eg.db",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

data_GO_sim <- simplify(data_GO, 
                        cutoff=0.7, 
                        by="p.adjust", 
                        select_fun=min)

dotplot(data_GO_sim, showCategory=5,font.size = 8)+
  theme(axis.text.x=element_text(angle=30,hjust=1))
data_GO_sim_fil <- data_GO_sim@compareClusterResult
dev.off()


###拟时序
table(scRNA_T.nmf$NMF_celltype)
library(BiocGenerics)
library(monocle)
library(tidyverse)
library(patchwork)
data=as.matrix(scRNA_T.nmf@assays$RNA@counts)

data <- as(data, 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNA_T.nmf@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
## 以下代码一律不得修???
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())

mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)

##使用monocle选择的高变基因
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
plot_ordering_genes(mycds)

#降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
mycds <- orderCells(mycds)
ss=intersect(gene,rownames(scRNA_T.nmf))
dev.off()
my_pseudotime_cluster <- plot_pseudotime_heatmap(mycds[ss,],
                                                 num_clusters = 6, # add_annotation_col = ac,
                                                 show_rownames = TRUE,
                                                 return_heatmap = TRUE)

my_pseudotime_cluster 
library(ClusterGVis)

plot_pseudotime_heatmap2(mycds[ss,],
                         num_clusters = 6,
                         cores = 1,
                         show_rownames = T,
                         return_heatmap = T)


df <- plot_pseudotime_heatmap2(mycds[ss,],
                               num_clusters = 6,
                               cores = 1)

visCluster(object = df,plot.type = "line")
visCluster(object = df,plot.type = "heatmap")

pdf(file = "results/拟时序.pdf",height = 8,width = 6,onefile = F)
gene1 = sample(df$wide.res$gene,30,replace = TRUE)
visCluster(object = df,plot.type = "heatmap",
           markGenes = gene1)

visCluster(object = df,plot.type = "both")
dev.off()

#State轨迹分布???
p1<- plot_cell_trajectory(mycds, color_by = "State")


p2 <- plot_cell_trajectory(mycds, color_by = "NMF_celltype")


p3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")

##合并出图
library(patchwork)#拼图
p2+p3
dev.off()

#接下来就是计算km曲线