
############################# Analysis of single cell RNA-seq data from Feldman et al. ###################################################

#Load requied libraries 
library(Seurat)
library(dplyr)
library(openxlsx)
library(UCell)
library(patchwork)
library(ggplot2)
library(reshape2)
library(reshape)
library(rstatix)
library(ggpubr)

# Load rds file Seurat object (Provided by TD and CK)
feldman <- readRDS("/Users/........../feldman_seurat.rds")

# Load annotation file 
ann <- read.xlsx("/Users/.........../Feldman_ann.xlsx", rowNames = TRUE)

# Add cell annotation to seurat object  
feldman <- AddMetaData(feldman, ann, col.name = NULL)

#NOT RUN # #when work on whole dataset: all cell types then for convenience dall dataset (feldman) is renamed as macro
#macro <- feldman


# Subset seurat object by taking the 
macro <- subset(feldman, subset = Cluster.number == "3")

# # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
macro[["percent.mt"]] <- PercentageFeatureSet(macro, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(macro, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(macro, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(macro, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Normalizing the data
#pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
macro <- NormalizeData(macro)

#Identification of highly variable features (feature selection)
macro <- FindVariableFeatures(macro, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(macro), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(macro)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Scaling the data
all.genes <- rownames(macro)
macro <- ScaleData(macro, features = all.genes)

#Perform linear dimensional reduction
macro <- RunPCA(macro, features = VariableFeatures(object = macro))


# Examine and visualize PCA results a few different ways
print(macro[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(macro, dims = 1:2, reduction = "pca")
macro_cluster <- RunUMAP(macro, dims = 1:10)
DimPlot(macro, reduction = "pca")

#PCA heatmap
DimHeatmap(macro, dims = 1, cells = 500, balanced = TRUE)

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
macro <- JackStraw(macro, num.replicate = 100)
macro <- ScoreJackStraw(macro, dims = 1:20)
JackStrawPlot(macro, dims = 1:15)


#Cluster the cells
macro <- FindNeighbors(macro, dims = 1:10)
macro <- FindClusters(macro, resolution = 0.5)

#UMAP only
macro <- RunUMAP(macro, dims = 1:10)
DimPlot(macro, reduction = "umap")

#UMAP with ann data cluster number
DimPlot(macro, reduction = "umap", group.by= "Cluster.number")

#save umap wth cluster no
ggsave("UMAP_Cluster.number2.pdf", width = 4, height = 4)

#Response UMAP
DimPlot(macro, reduction = "umap", group.by= "characteristics..response")

# Pre and Post
DimPlot(macro, reduction = "umap", group.by= "characteristics..patinet.ID..Pre.baseline..Post..on.treatment.")

#save rds
saveRDS(macro, file = "/Users/......../feldman_seurat_SC.rds")

#load saved RDS

#macro <- readRDS("feldman_seurat.rds")


# find all markers of cluster 0
cluster2.markers <- FindMarkers(macro, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 50)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones

macro.markers <- FindAllMarkers(macro, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
macro.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

#macro.markers <- pbmc.markers

#Save marker file for all clusters
write.xlsx(macro.markers, "/Users/......../Macro_clusters_Markers.xlsx")

#Show one marker in all clusters
VlnPlot(macro, features = c("TREM2",
                            "FRMD4A",
                            "PIK3IP1",
                            "CD209",
                            "RASGRP3",
                            "APOE"))

FeaturePlot(macro, features = c("TREM2",
                                "FRMD4A",
                                "PIK3IP1",
                                "CD209",
                                "RASGRP3", 
                                "APOE"))
#heatmap
macro.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(macro, features = top10$gene) + NoLegend()

# Meta data NR vs R
ICI_response <- macro$characteristics..response
ICI_response <- as.data.frame(ICI_response)

#################### END of INTIAL ANALYSIS ####################################

# load RDS file

macro <- readRDS("/Users/........../feldman_seurat_SC.rds")

#Load requied libraries 
library(Seurat)
library(dplyr)
library(openxlsx)
library(UCell)
library(patchwork)
library(ggplot2)

#Load gene set All TAM signature
TAM_sig = read.xlsx("/Users/........./LM22.xlsx", sheet = 4)
signatures <- lapply(1:ncol(TAM_sig), function(i) {as.character(TAM_sig[, i])} )
names(TAM_sig) <- colnames(TAM_sig)
TAM_sig_Ucell <- lapply(TAM_sig, function(i) unique(i[!is.na(i)]))


#Add Gene Signature to Seurat Onject 
macro <- AddModuleScore_UCell(macro, features = TAM_sig_Ucell)

#Save seurat.object with UCell score
saveRDS(macro, file = "feldman_seurat_UCell_SC.rds")

# saving seurat.object as data frame
ucll <- as.data.frame(macro[[]])
write.xlsx(ucll, '/Users/............/Feldmann_Ucell.xlsx')

#UCell map plots
signature.names <- paste0(names(TAM_sig_Ucell), "_UCell")

# UMAP
p <- FeaturePlot(macro, reduction = "umap", 
                 features = signature.names, ncol = 3, 
                 order = F, 
                 pt.size = 0.3,
                 cols = c("darkred", "orange", "yellow", "white"))
p

# UMAP for spefific signatures
p <- FeaturePlot(macro, reduction = "umap", 
                 features = c("IFN_TAMs_UCell", "LA_TAM_UCell"), ncol = 2,
                 split.by = "characteristics..response",
                 order = F, 
                 pt.size = 0.3,
                 cols = c("darkred", "orange", "yellow", "white"))
p

#Violin
v <- VlnPlot(macro, 
             features = signature.names,
             group.by = "seurat_clusters")
v

#load saved Ucell matrix
mat <- read.xlsx("/Users/.........../Tirosh_macro_Ucell.xlsx")


#clean mat
mat <- mat[, -4:-6]

# melt 
library(reshape2)
library(reshape)

mat_melt <- melt(mat)

#save long data

#write.xlsx(mat_melt, "Feldmann_Ucell_long_format.xlsx")


#violin zitter
library(ggplot2)

#Ucell score NR vs R
v <- ggplot(mat_melt, aes(x = variable, y = value, colour = variable)) +
  ylim(0, 0.4) +
  geom_jitter(position=position_jitter(0.2), size = 1, alpha = 0.5) +
  facet_wrap( variable ~ . ) +
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               colour = "black")+ 
  theme_bw()
v
#UCell in seurat clusters

v <- ggplot(mat_melt, aes(x = seurat_clusters, y = value, colour = seurat_clusters)) +
  ylim(0, 0.4) +
  geom_jitter(position=position_jitter(0.2), size = 1, alpha = 0.5) +
  facet_wrap( variable ~ . ) +
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               colour = "black")+ 
  theme_bw()
v
ggsave("Ucell_NR_R.pdf", width = 7, height = 5)




#stat
library(rstatix)
library(ggpubr)

stat.test <- mat_melt %>%
  group_by(variable) %>%
  t_test(value ~ characteristics..response) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test 

write.xlsx(stat.test, "ttest_Ucell_NR_R.xlsx")
ggsave("Ucell_NR_R.pdf", width = 7, height = 5)

################################################################################

#subseting datagrame based on pre and post

library(dplyr)
mat_cat <- split(mat, f = mat$characteristics..patinet.ID..Pre.baseline..Post..on.treatment.)                     # Split data
mat_cat                                                 # Print list

mat_cat_avg <- lapply(mat_cat, function(i){
  m <- i[, -c(1:7)]
  #return(apply(m, 2, median))
  return(colMeans(m))
})

mat_cat_avg <- do.call(rbind, mat_cat_avg)

#Covert to data grame 

cat_avg <- as.data.frame(mat_cat_avg)


# save avg ucell score across cell types

write.xlsx(cat_avg, "/Users/.........../Avg_UCell_CellTypes.xlsx", rowNames = TRUE)



# ggpaired
#data load

ggpair <- read.xlsx("/Users/........../Avg_UCell_CellTypes2_ggpaired_Format.xlsx")



# douple facet 
p <- ggpaired(ggpair, cond1 =  "Pre_LA_TAM_UCell", cond2 = "Post_LA_TAM_UCell",
              color = "condition", line.color = "gray", line.size = 0.4, alpha= 0.1,
              palette = "jco")+
  facet_wrap(~Response)
stat_compare_means(paired = TRUE)

p

write.xlsx(cat_avg, "/Users/......../Avg_UCell_CellTypes2.xlsx", rowNames = TRUE)

#stat
#load paired u cell matrix
paired_mat <- read.xlsx("/Users/......../Avg_UCell_CellTypes.xlsx")

dfNR <- paired_mat[paired_mat$Response == "NR", ]
dfR <- paired_mat[paired_mat$Response == "R", ]

#Data frame melt
dfNR_melt <- melt(dfNR)
dfR_melt <- melt(dfR)

dfNR_melt <- as.data.frame(dfNR_melt)
dfR_melt <- as.data.frame(dfR_melt)

#Stat
#NR
stat.test <- dfNR_melt %>%
  group_by(variable) %>%
  wilcox_test(value ~ Timing, paired = TRUE, alternative = "greater") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test 

write.xlsx(stat.test, "/Users/.........../NR_greater_Ucell_Wilcox.test.xlsx")

stat.test <- dfR_melt %>%
  group_by(variable) %>%
  t_test(value ~ Timing, paired = TRUE, alternative = "two.sided") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test 


write.xlsx(stat.test, "/Users/......../R_less_Ucell_Wilcox.test.xlsx")

###############################################################################

#Ucell score NR vs R and pre and post

melt_mat <- read.xlsx("/Users/.........../Feldmann_Ucell_long_format.xlsx")

v <- ggplot(melt_mat, aes(x = Timing, y = value, colour = Timing)) +
  ylim(0, 0.4) +
  geom_jitter(position=position_jitter(0.2), size = 1, alpha = 0.5) +
  facet_wrap( variable ~ Response) +
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               colour = "black")+ 
  theme_bw()
v1 <- v + scale_x_discrete(limits = c("Pre","Post"))

v1
ggsave("Pre_Post_NR_R_UCell_FINAL_PLOT.pdf", width = 8, height = 8)

#NR and R dataframe

dfNR <- melt_mat[melt_mat$Response == "Non-responder", ]
dfR <- melt_mat[melt_mat$Response == "Responder", ]

#Load melt library


#Data frame melt
dfNR_melt <- melt(dfNR)
dfR_melt <- melt(dfR)
dfNR_melt <- as.data.frame(dfNR_melt)
dfR_melt <- as.data.frame(dfR_melt)

#pre and post dataframe
df_pre <- melt_mat[melt_mat$Timing == "Pre", ]
df_post <- melt_mat[melt_mat$Timing == "Post", ]

dfNR_melt <- dfNR_melt[, -9]
dfR_melt <- dfR_melt[, -9]


#Stat
#NR: Pre vs. Post Greater
stat.test <- dfNR_melt %>%
  group_by(variable) %>%
  wilcox_test(value ~ Timing, paired = FALSE, alternative = "less") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test 

write.xlsx(stat.test, "/Users/............/NR_Pre_vs_Post_Less_Ucell_Wilcox.test.xlsx")

#R: Pre vs. Post Less
stat.test <- dfR_melt %>%
  group_by(variable) %>%
  wilcox_test(value ~ Timing, paired = FALSE, alternative = "less") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test 



#Stat
#NR
stat.test <- df_post %>%
  group_by(variable) %>%
  wilcox_test(value ~ Response, paired = FALSE, alternative = "greater") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test 

write.xlsx(stat.test, "/Users/......../Post_greater_Ucell_Wilcox.test.xlsx")

#R
stat.test <- dfR_melt %>%
  group_by(variable) %>%
  wilcox_test(value ~ Timing) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test 


write.xlsx(stat.test, "/Users/................./R_less_Ucell_Wilcox.test.xlsx")

#####################################################################################################################

# Preparing FILES for GSEA Clusterprofiler Input

Seurat_cluster_DEG <- read.xlsx("/Users/........./Macro_clusters_Markers.xlsx")

C0 <- Seurat_cluster_DEG[Seurat_cluster_DEG$cluster == 0, ]
C1 <- Seurat_cluster_DEG[Seurat_cluster_DEG$cluster == 1, ]  
C2 <- Seurat_cluster_DEG[Seurat_cluster_DEG$cluster == 2, ]
C3 <- Seurat_cluster_DEG[Seurat_cluster_DEG$cluster == 3, ]
C4 <- Seurat_cluster_DEG[Seurat_cluster_DEG$cluster == 4, ]
C5 <- Seurat_cluster_DEG[Seurat_cluster_DEG$cluster == 5, ]
C6 <- Seurat_cluster_DEG[Seurat_cluster_DEG$cluster == 6, ]
C7 <- Seurat_cluster_DEG[Seurat_cluster_DEG$cluster == 7, ]
C8 <- Seurat_cluster_DEG[Seurat_cluster_DEG$cluster == 8, ]

write.xlsx(C0, "C0.xlsx")
write.xlsx(C1, "C1.xlsx")
write.xlsx(C2, "C2.xlsx")
write.xlsx(C3, "C3.xlsx")
write.xlsx(C4, "C4.xlsx")
write.xlsx(C5, "C5.xlsx")
write.xlsx(C6, "C6.xlsx")
write.xlsx(C7, "C7.xlsx")
write.xlsx(C8, "C8.xlsx")

################################################################################
### Hyper G test on Seurat clusters ############################################
#Install hyperR
devtools::install_github("montilab/hypeR")
#load library
library(hypeR)

#From Seurat
macro.markers <- FindAllMarkers(macro, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
macro.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

#macro.markers <- pbmc.markers

#Save marker file for all clusters
write.xlsx(macro.markers, "/Users/........./Macro_clusters_Markers.xlsx")

# convert macro.markers dataframe as list 

list <- split(unlist(macro.markers$gene), macro.markers$cluster)
names(list) <- c("C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")

# msig DB all
C1 <- msigdb_gsets("Homo sapiens", "C1", clean=TRUE)
C2_KEGG <- msigdb_gsets("Homo sapiens", "C2", "CP:KEGG", clean=TRUE)
C2_Biocarta <- msigdb_gsets("Homo sapiens", "C2", "CP:BIOCARTA", clean=TRUE)
C2_REACTOME <- msigdb_gsets("Homo sapiens", "C2", "CP:REACTOME", clean=TRUE)
C2_WIKIPATHWAYS <- msigdb_gsets("Homo sapiens", "C2", "CP:WIKIPATHWAYS", clean=TRUE)
HALLMARK <- msigdb_gsets("Homo sapiens", "H", clean=TRUE)
C5_GOBP <- msigdb_gsets("Homo sapiens", "C5", "GO:BP", clean=TRUE)
C5_GOCC <- msigdb_gsets("Homo sapiens", "C5", "GO:CC", clean=TRUE)
C5_GOMF <- msigdb_gsets("Homo sapiens", "C5", "GO:MF", clean=TRUE)
C5_HPO <- msigdb_gsets("Homo sapiens", "C5", "HPO", clean=TRUE)
C6 <- msigdb_gsets("Homo sapiens", "C6", clean=TRUE)
C7 <- msigdb_gsets("Homo sapiens", "C7", "IMMUNESIGDB", clean=TRUE)
C8 <- msigdb_gsets("Homo sapiens", "C8", clean=TRUE)

# RUN HYPER G

C1_hyp <- hypeR(list, C1, test="hypergeometric", background=30000)
C2_KEGG_hyp <- hypeR(list, C2_KEGG, test="hypergeometric", background=30000)
C2_Biocarta_Hyp <- hypeR(list, C2_Biocarta, test="hypergeometric", background=30000)
C2_REACTOME_Hyp <- hypeR(list, C2_REACTOME, test="hypergeometric", background=30000)
C2_WIKI_Hyp <- hypeR(list, C2_WIKIPATHWAYS, test="hypergeometric", background=30000)
HALLMARK_hyp <- hypeR(list, HALLMARK, test="hypergeometric", background=30000)
GOBP_hyp <- hypeR(list, C5_GOBP, test="hypergeometric", background=30000)
GOMF_hyp <- hypeR(list, C5_GOMF, test="hypergeometric", background=30000)
GOCC_hyp <- hypeR(list, C5_GOCC, test="hypergeometric", background=30000)
HPO_hyp <- hypeR(list, C5_HPO, test="hypergeometric", background=30000)
C6_hyp <- hypeR(list, C6, test="hypergeometric", background=30000)
C7_hyp <- hypeR(list, C7, test="hypergeometric", background=30000)
C8_hyp <- hypeR(list, C8, test="hypergeometric", background=30000)

# plots
C1_plot <- hyp_dots(C1_hyp, merge=TRUE, fdr=0.05, title="C1")
C2_KEGG_plot <- hyp_dots(C2_KEGG_hyp, merge=TRUE, fdr=0.05, title="KEGG")
C2_BIOCARTA_plot <- hyp_dots(C2_Biocarta_Hyp, merge=TRUE, fdr=0.05, title="Biocarta")
C2_REACTOME_plot <- hyp_dots(C2_REACTOME_Hyp, merge=TRUE, fdr=0.05, title="REACTOME")
C2_wiki_plot <- hyp_dots(C2_WIKI_Hyp, merge=TRUE, fdr=0.05, title="WIKIPATHWAY")
Hallmark_plot <- hyp_dots(HALLMARK_hyp, merge=TRUE, fdr=0.05, title="HALLMARK")
GOBP_plot <- hyp_dots(GOBP_hyp, merge=TRUE, fdr=0.05, title="GOBP")
GOMF_plot <- hyp_dots(GOMF_hyp, merge=TRUE, fdr=0.05, title="GOMF")
GOCC_plot <- hyp_dots(GOCC_hyp, merge=TRUE, fdr=0.05, title="GOCC")
HPO_plot <- hyp_dots(HPO_hyp, merge=TRUE, fdr=0.05, title="HPO")
C6_plot <- hyp_dots(C6_hyp, merge=TRUE, fdr=0.05, title="C6")
C7_plot <- hyp_dots(C7_hyp, merge=TRUE, fdr=0.05, title="C7")
C8_plot <- hyp_dots(C8_hyp, merge=TRUE, fdr=0.05, title="C8")

#set dir
setwd("/Users/...../Melanoma/scRNA-seq/DEG/HyperG/")


#save plots 
ggsave(C1_plot, filename = "C1_plot.pdf",width = 8, height = 5)
ggsave(C2_KEGG_plot, filename = "KEGG_plot.pdf",width = 8, height = 5)
ggsave(C2_BIOCARTA_plot, filename = "BIOCARTA_plot.pdf",width = 8, height = 5)
ggsave(C2_REACTOME_plot, filename = "REACTOME.pdf",width = 8, height = 5)
ggsave(C2_wiki_plot, filename = "WIKIPATHWAYS.pdf",width = 8, height = 5)
ggsave(Hallmark_plot, filename = "HALLMARK.pdf",width = 8, height = 5)
ggsave(GOBP_plot, filename = "GOBP.pdf",width = 8, height = 5)
ggsave(GOMF_plot, filename = "GOMF.pdf",width = 8, height = 5)
ggsave(GOCC_plot, filename = "GOCC.pdf",width = 8, height = 5)
ggsave(HPO_plot, filename = "HPO.pdf",width = 8, height = 5)
ggsave(C6_plot, filename = "C6.pdf",width = 8, height = 5)
ggsave(C7_plot, filename = "C7.pdf",width = 8, height = 5)
ggsave(C8_hyp, filename = "C8.pdf",width = 8, height = 5)


############################  STAT on UCell scores across SEURAT CLUSTERS #######################

data <- read.xlsx("/Users/........./Feldmann_Ucell_long_format.xlsx")

#Stat

stat.test <- data %>%
  group_by(variable) %>%
  wilcox.test(value ~ seurat_clusters) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()
stat.test 

write.xlsx(stat.test, "/Users/........./Feldmann_Ucell_long_format_STAT.xlsx")
