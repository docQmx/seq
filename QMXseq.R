###################
#step1: create seurat object and doublet removal
#step2: data integration by harmony 
#step3: annotation by canonical markers

cat ("sample", i, "cell: ")

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())

library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(viridis)
library(DoubletFinder)
library(tidyverse)
library(devtools)
library(harmony)
library(MAST)
set.seed(1)

#################################################################################
#function definition
seurat.standard.normalize.and.scale <- function(seu_obj){
  seu_obj <- NormalizeData(seu_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 2000)
    setwd("F:/space/PLOT/")
    pdf(paste0("./VariableFeatures_plots_", project_name, ".pdf"))
    print(LabelPoints(plot = VariableFeaturePlot(seu_obj), points = head(VariableFeatures(seu_obj), 10), repel = TRUE, size=2.5) )
    dev.off()
  allgenes <- rownames(seu_obj)
  seu_obj <- ScaleData(seu_obj, features = allgenes)
  seu_obj <- RunPCA(seu_obj, features = VariableFeatures(object = seu_obj))
  seu_obj <- FindNeighbors(seu_obj, dims = 1:20)
  seu_obj <- FindClusters(seu_obj, resolution = 0.5)
  seu_obj <- RunUMAP(seu_obj, dims = 1:40)
  return(seu_obj)
}

make.seurat.object.and.doublet.removal <- function(data_directory, project_name){
  # function for basic seurat based qc and doubletfinder based doublet removal
  setwd("F:/space/data_mix/")
  data_counts <- Read10X(data.dir = data_directory)
  currentSample <- CreateSeuratObject(counts = data_counts, project = project_name, min.cells = 3, min.features = 40)
  currentSample[["percent.mt"]] <- PercentageFeatureSet(currentSample, pattern = "^MT-")
  # qc plot-pre filtering
  setwd("F:/space/PLOT/")
  pdf(paste0("./qc_plots_", project_name, "_prefiltered.pdf"))
  print(VlnPlot(currentSample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.05))
  dev.off()
  pdf(paste0("./qc_plots_", project_name, "_prefiltered_no_points.pdf"))
  print(VlnPlot(currentSample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0))
  dev.off()
  # select cells (nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt<5)
  currentSample <- subset(currentSample, subset =  nFeature_RNA > 300 & nCount_RNA > 1000 & percent.mt<5)
  # Normalize and make UMAP
  currentSample <- seurat.standard.normalize.and.scale(currentSample)
  # Run doublet finder (pN = 0.25 & PK calculation & 7.5% doublet formation rate)
  sweep.res.list_bc <- paramSweep(currentSample, PCs = 1:10, sct = FALSE)
  sweep.stats_bc <- summarizeSweep(sweep.res.list_bc, GT = FALSE)
  bcmvn_bc <- find.pK(sweep.stats_bc)
  mpK<-as.numeric(as.vector(bcmvn_bc$pK[which.max(bcmvn_bc$BCmetric)]))
  annotations <- currentSample@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)          
  nExp_poi <- round(0.075*nrow(currentSample@meta.data))
  seu_sample <- doubletFinder(currentSample, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = NULL, sct = FALSE)
    print(head(seu_sample@meta.data))
    seu_sample$doublet.class <- seu_sample[[paste0("DF.classifications_0.25_", mpK, "_", nExp_poi)]]
    seu_sample[[paste0("DF.classifications_0.25_", mpK, "_", nExp_poi)]] <- NULL
    pdf(paste0("./UAMP_doublet_", project_name, ".pdf"))
    print(DimPlot(seu_sample, reduction = "umap", group.by ="doublet.class",cols =c("black", "red")))
    dev.off()
  # UMAP plot-seurat clusters
  pdf(paste0("./UMAP_seuratclusters_", project_name, ".pdf"))
  print(DimPlot(seu_sample, reduction = "umap", group.by = "seurat_clusters" ))
  dev.off()
  # Remove extra stuff and return filtered Seurat object
  seu_sample <- DietSeurat(seu_sample, counts=TRUE, data=TRUE, scale.data=FALSE, assays="RNA")
    #setwd("F:/space/QMX/QMXseq/")
    #saveRDS(seu_sample,file = paste0(project_name, ".rds"))
  return(seu_sample)
}

seurat.qc.plots <- function(seu_obj, sample_name){
  # Make some basic qc plots
  pdf(paste0("./seurat_nFeature_plots_", sample_name, ".pdf"), width = 40, height = 15)
  print(VlnPlot(seu_obj, features = c("nFeature_RNA"), ncol = 1, pt.size = 0.2))
  dev.off()
  
  pdf(paste0("./seurat_nCount_plots_", sample_name, ".pdf"), width = 40, height = 15)
  print(VlnPlot(seu_obj, features = c("nCount_RNA"), ncol = 1, pt.size = 0.2))
  dev.off()
  
  pdf(paste0("./seurat_pMT_plots_", sample_name, ".pdf"), width = 40, height = 15)
  print(VlnPlot(seu_obj, features = c("percent.mt"), ncol = 1, pt.size = 0.2))
  dev.off()
}

#################################################################################
#data quality control and normalization
setwd("F:/space/data_mix/")
data_directory=c("Z","T3a","T3b","T2a","T4a","T4b","LN1","LN3a","LN3b","LN2a","LN4a2","LN4b2","LP2a","LP4a","LP4b")
project_name=c("T1","T2","T3","T4","T5","T6","LN1","LN2","LN3","LN4","LN5","LN6","LP4","LP5","LP6")
samples <- project_name
sample1 <- make.seurat.object.and.doublet.removal(data_directory[1], samples[1])
seu_list <- sample1
for (i in 2:length(samples)){
  sc.i = make.seurat.object.and.doublet.removal(data_directory[i], samples[i])
  seu_list=merge(seu_list,sc.i)
}
table(seu_list$orig.ident)
setwd("F:/space/QMX/QMXseq/")
saveRDS(seu_list,file = "seu_list.rds")

#data integration-harmony
scRNA_harmony=seu_list
scRNA_harmony  <- NormalizeData(scRNA_harmony ) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution =0.5)
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:20)
scRNA_harmony$orig.ident <- factor(scRNA_harmony$orig.ident,levels = project_name)
saveRDS(scRNA_harmony,file = "scRNA_harmony.rds")
table(scRNA_harmony@meta.data$seurat_clusters)
DimPlot(scRNA_harmony, reduction = "umap",label = T) 
DimPlot(scRNA_harmony, reduction = "umap", split.by ='orig.ident',ncol = 5)
DimPlot(scRNA_harmony, reduction = "umap", group.by='orig.ident')
table(scRNA_harmony1$orig.ident)  

#annotation
##findermarker
Idents(scRNA_harmony)="seurat_clusters"
table(scRNA_harmony1@meta.data$seurat_clusters)
scRNA_anno=scRNA_harmony
scRNA_anno<- JoinLayers(scRNA_anno)
table(scRNA_anno@active.ident)
markers <- FindAllMarkers(scRNA_anno, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 1,test.use = "MAST")
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
##canonical markers
###XUNYIN
T_cells = c("CD247","THEMIS", "SKAP1", "CD96","ITK")
B_cells = c("CD79A", "CD79B", "MS4A1")
Plasma_cells = c("JCHAIN","IGHG1","IGKC","IGLC2")
Myeloid_cells= c("C1QB","C1QC","CAP3","CD14","CD163")
DCs = C("LAMP3","CLEC9A")
Fibroblasts = c("COL14A1", "COL1A1", "COL6A3","COL5A2")
Endothelial_cells = c("PCAM1", "VWF")
Epithelial_cells = c("EPCAM", "KRT18","KRT8","KRT14","KRT5","MID1","ANKS6")
NK_cells = c("NCR1","GLNY","NKG7") 
###Nature
T_cells = c("ALOX5AP","BST2","CCL3","CD4",
            "CD40LG","CD82","CD8A","CXCL13","DUSP1","ENTPD1",
            "GNLY","GPR25","GZMB","GZMH","HAVCR2","HLA-DQA2",
            "HLA-DRA","HLA-DRB1","HLA-DRB5","IFI27","IFI44",
            "IFI44L","IFI6","IFIT3","IFITM3","IL17","IL7R",
            "INPP5F","ISG15","KRT86","LY6E","MX1","OAS1","PRF1")
B_cells = c("BANK1","CCDC50","CD14","CD19","CD27","CD38",
            "CD40","CD74","CD79A","IGHD","IGHM","IGJ",
            "IGKC","IGLL1","IGLL5","LAPTM5","MS4A1","PTPRC",
            "SEL1L3","SSR4","TCF4","TMP1","VIM")
###TOP
#B_cells = c("CD79A","MS4A1","BANK1","BLK")
#T_cells = c("CD247","THEMIS", "CD96","ITK")
#Plasma_cells = c("JCHAIN","IGHG1","IGKC","IGLC2")
#NK_cells = c("NCR1","GLNY","NKG7")
Lymphoid_cells=C("CD79A","MS4A1","BANK1","BLK","CD247","THEMIS", "CD96","ITK")
Fibroblasts = c("COL14A1", "COL1A1", "COL6A3","COL5A2")
Endothelial_cells = c("VWF","FLT1","CDH5","PLVAP")
Epithelial_cells = c("EPCAM","KRT18","KRT8","ANKRD30A")
Neurons=c("NRXN1","SNAP25","MAP2","KIF5A")
Monocytes=c("HLA-DRA","HLA-DRB1","CD300E","FCN1")
Macrophages=c("CD163","C1QB","C1QC","APOE")
Dendriti_cells= c("LAMP3","CCR7","FSCN1","LY75")

Idents(scRNA_harmony)="seurat_clusters"
DotPlot(scRNA_harmony,features =c("EPCAM","KRT18","KRT8","ANKRD30A"))+RotatedAxis()
FeaturePlot(scRNA_harmony,features =c("EPCAM","KRT18","KRT8","ANKRD30A"))
VlnPlot(scRNA_harmony,
        features = c("EPCAM","KRT18","KRT8","ANKRD30A"),
        group.by = "seurat_clusters",
        pt.size = 0, 
        ncol = 2)
#scRNA_filtered <- subset(scRNA_harmony, subset = !is.na(CD300E) & !is.na(FCN1))
#VlnPlot(scRNA_filtered,features = c("HLA-DRA","HLA-DRB1","CD300E","FCN1"),group.by = "seurat_clusters",pt.size = 0, ncol = 2)

scRNA_renameharmony=RenameIdents(scRNA_harmony,
                            "0"="Epithelial Cells",
                            "1"="Epithelial Cells" ,
                            "2"="Fibroblasts",
                            "3"="Epithelial Cells",
                            "4"="Lymphoid Cells", 
                            "5"="Macrophages" ,
                            "6"="Endothelial Cells",
                            "7"="Epithelial Cells",
                            "8"="Lymphoid Cells",
                            "9"="Endothelial Cells",
                           "10"="Monocytes",
                           "11"="Epithelial Cells",
                           "12"="Dendritic Cells",
                           "13"="Epithelial Cells",
                           "14"="Lymphoid Cells",
                           "15"="Epithelial Cells",
                           "16"="Neurons",
                           "17"="Epithelial Cells")
DimPlot(scRNA_renameharmony,label = T)
scRNA_harmony@meta.data$celltype=scRNA_renameharmony@active.ident
colnames(scRNA_renameharmony@meta.data)
DimPlot(scRNA_renameharmony,group.by ="celltype",split.by = "orig.ident",ncol = 6)

Idents(scRNA_harmony)="orig.ident"
Idents(scRNA_harmony)="seurat_clusters"
genes_to_check <- list(Epithelial_cells = c("EPCAM", "KRT18", "KRT8", "ANKRD30A"),
                       Fibroblasts = c("COL14A1", "COL1A1", "COL6A3", "COL5A2"),
                       Lymphoid_cells = c("CD79A", "MS4A1", "BANK1", "BLK", 
                                          "CD247", "THEMIS", "CD96", "ITK"),
                       Macrophages = c("CD163", "C1QB", "C1QC", "APOE"),
                       Endothelial_cells = c("VWF", "FLT1", "CDH5", "PLVAP"),
                       Monocytes = c("HLA-DRA", "HLA-DRB1", "CD300E", "FCN1"),
                       Dendriti_cells = c("LAMP3", "CCR7", "FSCN1", "LY75"),
                       Neurons = c("NRXN1", "SNAP25", "MAP2", "KIF5A"))
DotPlot(scRNA_renameharmony, features = genes_to_check) +
  theme(axis.title = element_blank(),strip.text.x = element_blank(),axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.5))

features <- c("EPCAM", "KRT18", "KRT8", "ANKRD30A",
              "COL14A1", "COL1A1", "COL6A3", "COL5A2",
              "CD79A", "MS4A1", "BANK1", "BLK", 
              "CD247", "THEMIS", "CD96", "ITK",
              "CD163", "C1QB", "C1QC", "APOE",
              "VWF", "FLT1", "CDH5", "PLVAP",
              "HLA-DRA", "HLA-DRB1", "CD300E", "FCN1",
              "LAMP3", "CCR7", "FSCN1", "LY75",
              "NRXN1", "SNAP25", "MAP2", "KIF5A")
features <- c("EPCAM", "KRT18","KRT8","KRT14","KRT5","MID1","ANKS6")
scRNA_harmony=ScaleData(scRNA_harmony,features = rownames(scRNA_harmony))
DoHeatmap( scRNA_harmony, features = features, size = 3)



