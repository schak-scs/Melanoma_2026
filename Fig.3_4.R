# # install spata 2.0
# install.packages("devtools")
# library(BiocManager)
# 
# if (!base::requireNamespace("BiocManager", quietly = TRUE)){
#   install.packages("BiocManager")
# }
# 
# BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
#                        'limma', 'S4Vectors', 'SingleCellExperiment',
#                        'SummarizedExperiment', 'batchelor', 'Matrix.utils', 'EBImage'))
# 
# install.packages("Seurat")
# 
# # install tensorflow
# devtools::install_github(repo = "kueckelj/confuns")
# devtools::install_github(repo = "theMILOlab/SPATAData")
# devtools::install_github(repo = "theMILOlab/SPATA2")
# 
# # install.packages("remotes")
# remotes::install_github("rstudio/tensorflow")
# reticulate::install_python()
# 
# library(tensorflow)
# install_tensorflow(envname = "r-tensorflow")
# 
# install.packages("keras")
# library(keras)
# install_keras()
# 
# library(tensorflow)

# # if you want to use monocle3 related wrappers 
# devtools::install_github('cole-trapnell-lab/leidenbase')
# devtools::install_github('cole-trapnell-lab/monocle3')


# load library
library(SPATA2)
library(openxlsx)
library(ggplot2)
library(devtools)
#library(monocle3)
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(infercnv)
library(Seurat)
#library(SeuratData)
library(patchwork)
library(dplyr)
#library(dplyr)


# install BiocManager::install('glmGamPoi')
#BiocManager::install("glmGamPoi")

############################ SEURAT OBJECT ####################################

#seurat from 10x
data_dir = "/S1"
list.files(data_dir) # Should show filtered_feature_bc_matrix.h5
Seurat_S1 <- Load10X_Spatial(data.dir = data_dir)

# data processing 
plot1 <- VlnPlot(Seurat_S1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(Seurat_S1, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

################### Crteate a "data" slot in seurat object

seurat_object <- Load10X_Spatial(data.dir = "S1")

object <- seurat_object
object[["RNA3"]] <- as(object = object[["Spatial"]], Class = "Assay")
DefaultAssay(object) <- "RNA3"
object[["Spatial"]] <- NULL

object <- RenameAssays(object = object, RNA3 = 'Spatial')

spata_object <- SPATA2::asSPATA2(
  object = object,
  assay_name = "Spatial",
  sample_name = "S1",
  image_name = "slice1", 
  spatial_method = "Visium"
)

#convert seurat to spata2
spata_obj <- transformSeuratToSpata(object,
                                    sample_name = "Mel", method = "spatial", 
                                    assay_name = "Spatial")

# laod spata object 

S1 <- loadSpataObject("S1.rds")
S2 <- loadSpataObject("S2.rds")
S3 <- loadSpataObject("S3.rds")
S4 <- loadSpataObject("S4.rds")

plotSurfaceInteractive(object = S1)

#Autoencoder Denoising

#S1
# all expression matrices before denoising
getExpressionMatrixNames(object = S1)
# active expression matrix before denoising
getActiveMatrixName(object = S1)

# denoising your data 
S1 <-
  runAutoencoderDenoising(
    object = S1, 
    activation = "selu", 
    bottleneck = 56, 
    epochs = 20, 
    layers = c(128, 64, 32), 
    dropout = 0.1
  )

# all expression matrices after denoising
getExpressionMatrixNames(object = S1)

#S2
# all expression matrices before denoising
getExpressionMatrixNames(object = S2)
# active expression matrix before denoising
getActiveMatrixName(object = S2)

# denoising your data 
S2 <-
  runAutoencoderDenoising(
    object = S2, 
    activation = "selu", 
    bottleneck = 56, 
    epochs = 20, 
    layers = c(128, 64, 32), 
    dropout = 0.1
  )

# all expression matrices after denoising
getExpressionMatrixNames(object = S2)

#S3
# all expression matrices before denoising
getExpressionMatrixNames(object = S3)
# active expression matrix before denoising
getActiveMatrixName(object = S3)

# denoising your data 
S3 <-
  runAutoencoderDenoising(
    object = S3, 
    activation = "selu", 
    bottleneck = 56, 
    epochs = 20, 
    layers = c(128, 64, 32), 
    dropout = 0.1
  )

# all expression matrices after denoising
getExpressionMatrixNames(object = S3)

#S4
# all expression matrices before denoising
getExpressionMatrixNames(object = S4)
# active expression matrix before denoising
getActiveMatrixName(object = S4)

# denoising your data 
S4 <-
  runAutoencoderDenoising(
    object = S4, 
    activation = "selu", 
    bottleneck = 56, 
    epochs = 20, 
    layers = c(128, 64, 32), 
    dropout = 0.1
  )


#S4
# all expression matrices before denoising
getExpressionMatrixNames(object = spata_object)
# active expression matrix before denoising
getActiveMatrixName(object = spata_object)

# denoising your data 
spata_object <-
  runAutoencoderDenoising(
    object = spata_object, 
    activation = "selu", 
    bottleneck = 56, 
    epochs = 20, 
    layers = c(128, 64, 32), 
    dropout = 0.1
  )


# all expression matrices after denoising
getExpressionMatrixNames(object = S4)
getExpressionMatrixNames(object = spata_object)

############################# OBJECT CTEATIONA, DENOISING AND SAVE COMPLETE ###############################################

############################### CLUSTERING ################################################################################

# S1 
plotImageGgplot(object = S1)

# Check current clustering 
# current grouping options
getGroupingOptions(S1)

# run new clustering

# run the pipeline
S1 <- 
  runBayesSpaceClustering(
    object = S1, 
    name = "bayes_space" # the name of the output grouping variable
  )

# results are immediately stored in the objects feature data
getGroupingOptions(S1)

# plot
plotSurface(
  object = S1, 
  color_by = "bayes_space", 
  pt_clrp = "uc",
  display_image = FALSE
)

# S2 
plotImageGgplot(object = S2)

# Check current clustering 
# current grouping options
getGroupingOptions(S2)

# run new clustering
# run the pipeline
S2 <- 
  runBayesSpaceClustering(
    object = S2, 
    name = "bayes_space", # the name of the output grouping variable
    )


# results are immediately stored in the objects feature data
getGroupingOptions(S2)

# plot
plotSurface(
  object = S2, 
  color_by = "bayes_space", 
  pt_clrp = "uc",
  display_image = FALSE
)

# S3 
plotImageGgplot(object = S3)

# Check current clustering 
# current grouping options
getGroupingOptions(S3)

# run new clustering
# run the pipeline
S3 <- 
  runBayesSpaceClustering(
    object = S3, 
    name = "bayes_space" # the name of the output grouping variable
  )

# results are immediately stored in the objects feature data
getGroupingOptions(S3)

# plot
plotSurface(
  object = S3, 
  color_by = "bayes_space", 
  pt_clrp = "uc",
  display_image = FALSE
)

# S4 
plotImageGgplot(object = S4)

# Check current clustering 
# current grouping options
getGroupingOptions(S4)

# run new clustering
# run the pipeline
S4 <- 
  runBayesSpaceClustering(
    object = S4, 
    name = "bayes_space", # the name of the output grouping variable
    q_force = 5, # if need be desired number of cluster can be imposed 
    overwrite = TRUE
  )

# results are immediately stored in the objects feature data
getGroupingOptions(S4)

# plot
plotSurface(
  object = S4, 
  color_by = "bayes_space", 
  pt_clrp = "uc",
  display_image = FALSE
)

########################### END OF CLUSTERING ######################################

################################## DEA across clusters #############################################################

# check grouping option
getGroupingOptions(object = S1)

#Running the analysis
S1 <- runDEA(object = S1, across = "bayes_space", method_de = "wilcox")
S2 <- runDEA(object = S2, across = "bayes_space", method_de = "wilcox")
S3 <- runDEA(object = S3, across = "bayes_space", method_de = "wilcox")
S4 <- runDEA(object = S4, across = "bayes_space", method_de = "wilcox")

#check results 
printDeaOverview(object = S1)
printDeaOverview(object = S2)
printDeaOverview(object = S3)
printDeaOverview(object = S4)

# Extracting results

# extract the complete data.frame
#S1
S1_DEA <- getDeaResultsDf(
          object = S1, 
          across = "bayes_space", 
          method_de = "wilcox",
          n_highest_lfc = 500, # top 500 genes
          max_adj_pval = 0.01
          )
#S2
S2_DEA <- getDeaResultsDf(
  object = S2, 
  across = "bayes_space", 
  method_de = "wilcox",
  n_highest_lfc = 500, # top 500 genes
  max_adj_pval = 0.01
  )

write.xlsx(S2_DEA, "/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATA_object/S2/Cluster_DEA_ligand_receptor/S2_DEA.xlsx")

#S3
S3_DEA <- getDeaResultsDf(
  object = S3, 
  across = "bayes_space", 
  method_de = "wilcox",
  n_highest_lfc = 500, # top 500 genes
  max_adj_pval = 0.01
)



#S4
S4_DEA <- getDeaResultsDf(
  object = S4, 
  across = "bayes_space", 
  method_de = "wilcox",
  n_highest_lfc = 500, # top 500 genes
  max_adj_pval = 0.01
)

######################################  VISALIZATION OF DEA ##################################################################

# Heatmap

#S1
hmS1 <- 
  plotDeaHeatmap(
    object = S1, 
    across = "bayes_space",
    method_de = "wilcox",
    n_highest_lfc = 10, 
    n_bcsp = 100
  )

hmS1

#S2
hmS2 <- 
  plotDeaHeatmap(
    object = S2, 
    across = "bayes_space",
    method_de = "wilcox",
    n_highest_lfc = 10, 
    n_bcsp = 100
  )

hmS2

#S3
hmS3 <- 
  plotDeaHeatmap(
    object = S3, 
    across = "bayes_space",
    method_de = "wilcox",
    n_highest_lfc = 10, 
    n_bcsp = 100
  )

hmS3

#S4
hmS4 <- 
  plotDeaHeatmap(
    object = S4, 
    across = "bayes_space",
    method_de = "wilcox",
    n_highest_lfc = 10, 
    n_bcsp = 100
  )

hmS4


############################### END OF DEA #############################################################

############################### START OF LIGAND RECEPTOR DATABASE MAPPING ##############################
#INSTALL omnipath
install_github('saezlab/OmnipathR')

# laod
library(OmnipathR)
## We check some of the different interaction databases
get_interaction_resources()
## We check some of the different intercell categories
get_intercell_generic_categories()
## We import the intercell data into a dataframe
intercell <- import_omnipath_intercell(
  scope = "generic",
  aspect = "locational"
)

## We check the intercell annotations for the individual components of
## our previous complex. We filter our data to print it in a good format
S2_DEA_mapping <- dplyr::filter(intercell, genesymbol %in% S2_DEA$gene) %>%
  dplyr::distinct(genesymbol, parent, .keep_all = TRUE) %>%
  dplyr::select(category, genesymbol, parent) %>%
  dplyr::arrange(genesymbol)

#high conf 
icn <- import_intercell_network(high_confidence = TRUE)

## We check the intercell annotations for the individual components of
## our previous complex. We filter our data to print it in a good format
S2_DEA_mapping <- dplyr::filter(icn, icn$source_genesymbol %in% S2_DEA$gene) %>%
  dplyr::distinct(icn$source_genesymbol, parent, .keep_all = TRUE) %>%
  dplyr::arrange(icn$source_genesymbol)

############################### START OF HYPERG ########################################################

#load library
library(hypeR)
library(openxlsx)


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

# now define list_cluster_markers
#list_cluster_markers <- list_S1
#list_cluster_markers <- list_S2
#list_cluster_markers <- list_S3
list_cluster_markers <- list_S4

# RUN HYPER G for cluster markers

C1_hyp_cluster <- hypeR(list_cluster_markers, C1, test="hypergeometric", background=30000)
C2_KEGG_hyp_cluster <- hypeR(list_cluster_markers, C2_KEGG, test="hypergeometric", background=30000)
C2_Biocarta_Hyp_cluster <- hypeR(list_cluster_markers, C2_Biocarta, test="hypergeometric", background=30000)
C2_REACTOME_Hyp_cluster <- hypeR(list_cluster_markers, C2_REACTOME, test="hypergeometric", background=30000)
C2_WIKI_Hyp_cluster <- hypeR(list_cluster_markers, C2_WIKIPATHWAYS, test="hypergeometric", background=30000)
HALLMARK_hyp_cluster <- hypeR(list_cluster_markers, HALLMARK, test="hypergeometric", background=30000)
GOBP_hyp_cluster <- hypeR(list_cluster_markers, C5_GOBP, test="hypergeometric", background=30000)
GOMF_hyp_cluster <- hypeR(list_cluster_markers, C5_GOMF, test="hypergeometric", background=30000)
GOCC_hyp_cluster <- hypeR(list_cluster_markers, C5_GOCC, test="hypergeometric", background=30000)
HPO_hyp_cluster <- hypeR(list_cluster_markers, C5_HPO, test="hypergeometric", background=30000)
C6_hyp_cluster <- hypeR(list_cluster_markers, C6, test="hypergeometric", background=30000)
C7_hyp_cluster <- hypeR(list_cluster_markers, C7, test="hypergeometric", background=30000)
C8_hyp_cluster <- hypeR(list_cluster_markers, C8, test="hypergeometric", background=30000)


# plots for cell and clsuter # change directory to save different plots

C1_plot <- hyp_dots(C1_hyp_cluster, merge=TRUE, pval=0.05, title="C1")
C2_KEGG_plot <- hyp_dots(C2_KEGG_hyp_cluster, merge=TRUE, pval=0.05, title="KEGG")
C2_BIOCARTA_plot <- hyp_dots(C2_Biocarta_Hyp_cluster, merge=TRUE, pval=0.05, title="Biocarta")
C2_REACTOME_plot <- hyp_dots(C2_REACTOME_Hyp_cluster, merge=TRUE, pval=0.05, title="REACTOME")
C2_wiki_plot <- hyp_dots(C2_WIKI_Hyp_cluster, merge=TRUE, pval=0.05, title="WIKIPATHWAY")
Hallmark_plot <- hyp_dots(HALLMARK_hyp_cluster, merge=TRUE, pval=0.05, title="HALLMARK")
GOBP_plot <- hyp_dots(GOBP_hyp_cluster, merge=TRUE, pval=0.05, title="GOBP")
GOMF_plot <- hyp_dots(GOMF_hyp_cluster, merge=TRUE, pval=0.05, title="GOMF")
GOCC_plot <- hyp_dots(GOCC_hyp_cluster, merge=TRUE, pval=0.05, title="GOCC")
HPO_plot <- hyp_dots(HPO_hyp_cluster, merge=TRUE, pval=0.05, title="HPO")
C6_plot <- hyp_dots(C6_hyp_cluster, merge=TRUE, pval=0.05, title="C6")
C7_plot <- hyp_dots(C7_hyp_cluster, merge=TRUE, pval=0.05, title="C7")
C8_plot <- hyp_dots(C8_hyp_cluster, merge=TRUE, pval=0.05, title="C8")


#set dir


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
ggsave(C8_plot, filename = "C8.pdf",width = 8, height = 5)

#Save results
hyp_to_excel(C1_hyp_cluster, file_path="C1.xlsx")
hyp_to_excel(C2_KEGG_hyp_cluster, file_path="KEGG.xlsx")
hyp_to_excel(C2_Biocarta_Hyp_cluster, file_path="BIOCARTA.xlsx")
hyp_to_excel(C2_REACTOME_Hyp_cluster, file_path="REACTOME.xlsx")
hyp_to_excel(C2_WIKI_Hyp_cluster, file_path="WIKI.xlsx")
hyp_to_excel(HALLMARK_hyp_cluster, file_path="HALLMARK.xlsx")
hyp_to_excel(GOBP_hyp_cluster, file_path="GOBP.xlsx")
hyp_to_excel(GOMF_hyp_cluster, file_path="GOMF.xlsx")
hyp_to_excel(GOCC_hyp_cluster, file_path="GOCC.xlsx")
hyp_to_excel(HPO_hyp_cluster, file_path="HPO.xlsx")
hyp_to_excel(C6_hyp_cluster, file_path="C6.xlsx")
hyp_to_excel(C7_hyp_cluster, file_path="C7.xlsx")
hyp_to_excel(C8_hyp_cluster, file_path="C8.xlsx")

###############################################################################################################

######################################  GENE SET ENRICHMENT AND SURFACE PLOT ##################################

###############################################################################################################

# set spata object 
#spata_obj <- S1
#spata_obj <- S2 
#spata_obj <- S3
spata_obj <- S4


# GSEA
# DEFINE TAM SIGNATURES
geneMat_TAM <- read.xlsx(GeneSet.xlsx"))
geneList_TAM <- lapply(1:ncol(geneMat_TAM), function(i) as.character(geneMat_TAM[, i]))
names(geneList_TAM) <- colnames(geneMat_TAM)
geneList_TAM <- lapply(geneList_TAM, function(i) unique(i[!is.na(i)]))
# geneList$ALL <- unique(unlist(geneList))

# # DEFINE Melanoma SIGNATURES
# geneMat <- read.xlsx(file.path("Tsoi_2018_Differentiation.xlsx"), sheet = 1)
# geneList <- lapply(1:ncol(geneMat), function(i) as.character(geneMat[, i]))
# names(geneList) <- colnames(geneMat)
# geneList <- lapply(geneList, function(i) unique(i[!is.na(i)]))
# # geneList$ALL <- unique(unlist(geneList))



# ############################
# # Basic extracting functions
# 
# if(FALSE){
#   # the essential data.frame
#   getSpataDf(spata_obj)
#   
#   # dimensional reduction data
#   getUmapDf(spata_obj)
#   
#   # barcode spot coordinates
#   getCoordsDf(spata_obj)
# }
# 
# coords_df <- getCoordsDf(spata_obj)
# coords_df
# 
# joinWith(object = spata_obj, 
#          spata_df = coords_df,
#          features = "seurat_clusters", # cluster belonging
#          verbose = FALSE)
# 
# # output 
# joined_df


# GENES
#spata.genes.meta <- getGeneMetaData(spata_obj)
spata.genes <- getGenes(spata_obj)
geneList_TAM <- lapply(geneList_TAM, intersect, y = spata.genes)


### Adding TAM genesets
spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'IFN_TAMs', genes = geneList_TAM$IFN_TAMs, overwrite = TRUE)
#spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'Inflam_TAM', genes = geneList_TAM$Inflam_TAM)
spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'LA_TAM', genes = geneList_TAM$LA_TAM, overwrite = TRUE)
# spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'Angio_TAM', genes = geneList_TAM$Angio_TAM)
# spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'Reg_TAM', genes = geneList_TAM$Reg_TAM)
# spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'Prolif_TAM', genes = geneList_TAM$Prolif_TAM)
# spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'RTM_TAM', genes = geneList_TAM$RTM_TAM)
# spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'Macrophages.M2', genes = geneList_TAM$Macrophages.M2)
# spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'Macrophages.M1', genes = geneList_TAM$Macrophages.M1)
#spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'Undifferentiated-Neural.crest-like', genes = geneList_TAM$`Undifferentiated-Neural.crest-like`)
spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'Neural.crest-like', genes = geneList_TAM$`Neural.crest-like`, overwrite = TRUE)
spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'Transitory', genes = geneList_TAM$Transitory, overwrite = TRUE)
spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'Transitory-Melanocytic', genes = geneList_TAM$`Transitory-Melanocytic`, overwrite = TRUE)
spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'Melanocytic', genes = geneList_TAM$Melanocytic, overwrite = TRUE)
spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'Neural.crest-like-Transitory', genes = geneList_TAM$`Neural.crest-like-Transitory`, overwrite = TRUE)
spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'Exhausion_cell_cycle_CD8', genes = geneList_TAM$Exhausion_cell_cyle_CD8, overwrite = TRUE)
spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'Exhausted_CD8', genes = geneList_TAM$Exhausted_CD8, overwrite = TRUE)
#spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'Exhausion_HSP', genes = geneList_TAM$Exhausion_HSP)
spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'Memory_effector_1', genes = geneList_TAM$Memory_effector_1, overwrite = TRUE)
spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'Memory_effector_2', genes = geneList_TAM$Memory_effector_2, overwrite = TRUE)
spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'Activated_CD8', genes = geneList_TAM$Activated_CD8, overwrite = TRUE)
spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'Treg', genes = geneList_TAM$Treg, overwrite = TRUE)
#spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'Activated_CD8', genes = geneList_TAM$Exhausion_cell_cyle_CD8, overwrite = TRUE)
# spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'AXL_Melnoma', genes = geneList_TAM$AXL_Melnoma)
# spata_obj <- addGeneSet(object = spata_obj, class_name = 'mygs', gs_name = 'MITF_Melanoma', genes = geneList_TAM$MITF_Melanoma)

# save spata objects with gene sets

# set/change the current default directory
spata_obj <- setSpataDir(spata_obj, dir = "/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATA_object/S1/S1_denoised_GeneSet.RDS")

# save all denoised spata object 
saveSpataObject(object = spata_obj)

# set/change the current default directory
spata_obj <- setSpataDir(spata_obj, dir = "/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATA_object/S2/S2_denoised_GeneSet.RDS")

# save all denoised spata object 
saveSpataObject(object = spata_obj)

# set/change the current default directory
spata_obj <- setSpataDir(spata_obj, dir = "/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATA_object/S3/S3_denoised_GeneSet.RDS")

# save all denoised spata object 
saveSpataObject(object = spata_obj)

# set/change the current default directory
spata_obj <- setSpataDir(spata_obj, dir = "/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATA_object/S4/S4_denoised_GeneSet.RDS")

# save all denoised spata object 
saveSpataObject(object = spata_obj)



##############
# PLOT SURFACE


#plots$S4_Image

genesetS2 <- c("mygs_IFN_TAMs",
               "mygs_Melanocytic",
               "mygs_Memory_effector_1",
               "mygs_LA_TAM",
               "mygs_Exhausted_CD8",
               "mygs_Neural.crest-like",
               "mygs_Transitory-Melanocytic",
               "mygs_Transitory",
               "mygs_Activated_CD8",
               "mygs_Memory_effector_2",
               "mygs_Treg"
              )


# check which is the active matrix
getActiveExpressionMatrixName(spata_obj)

# if needed set the active matrix
spata_obj <- setActiveExpressionMatrix(spata_obj, "denoised")
spata_obj <- setActiveExpressionMatrix(spata_obj, "scaled")


# surface plot 

# open application to obtain a list of plots
#plots <- plotSurfaceInteractive(object = spata_obj)

# ssGSEA
S4_GS <- plotSurfaceComparison(object = spata_obj, 
                           color_by = genesetS2,
                           method_gs = "ssgsea",
                           smooth = TRUE,
                           pt_clrsp = "inferno",
                           smooth_span = 0.2,
                           pt_size = 2,
                           display_image = FALSE
                           )

S4_GS

ggsave(plot = S3_GS, file=file.path("/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATA_object/S4/Gene_Set_Enrichment/zscore/S4_zscore_genes_Surface_plot_newGS.pdf"),
       width = 6, height = 5)


GOI <- c("APOE", "TREM2", "IFNG", "IL10", "IL10RA", "IL10RB", "LDLR", "VLDLR", "LRP1", "LRP5", "LRP8", "AXL", "MITF", "TGFB1", "TGFBR1", "TGFBR2")

TGF <- c("APP",
         "ATF2",
         "CITED1",
         "ENG",
         "FOS",
         "FOSB",
         "ITGAV",
         "ITGB3",
         "JUN",
         "JUNB",
         "KLF10",
         "LEF1",
         "MMP1",
         "NUP153",
         "PML",
         "PPP1CB",
         "SPTBN1",
         "STAT1",
         "STAT3",
         "TGFB1",
         "TGFBR1",
         "TGFBR2",
         "THBS1",
         "TNC",
         "UBB",
         "UBE2I")



IFN <- c("CASP1",
         "CASP4",
         "CCL2",
         "CCL3",
         "CCL4",
         "CCL7",
         "CCL8",
         "CD274 ",
         "CD40",
         "CXCL2",
         "CXCL3",
         "CXCL9",
         "CXCL10",
         "CXCL11",
         "IDO1",
         "IFI6",
         "IFIT1",
         "IFIT2",
         "IFIT3",
         "IFITM1",
         "IFITM3",
         " IRF1",
         "IRF7",
         "ISG15",
         "LAMP3",
         "PDCD1LG2",
         "TNFSF10",
         "C1QA",
         "CD38",
         "IL4I1",
         "ISG15",
         "TNFSF10",
         "IFI44L")

LATAM <- c("ACP5",
           "APOE",
           "APOC1",
           "ATF1",
           "C1QA",
           "C1QB",
           "C1QC",
           "CCL18",
           "CD163",
           "CD36",
           "CD63",
           "CHI3L1",
           "CTSB",
           "CTSD",
           "CTSL",
           "F13A1",
           "FABP5",
           "FOLR2",
           "GPNMB",
           "IRF3",
           "LGALS3",
           "LIPA",
           "LPL",
           "MACRO",
           "MerTK",
           "MMP7",
           "MMP9",
           "MMP12",
           "MRC1",
           "NR1H3",
           "NRF1",
           "NUPR1",
           "PLA2G7",
           "RNASE1",
           "SPARC",
           "SPP1",
           "TFDP2",
           "TREM2",
           "ZEB1")


S1_exh <- plotSurfaceComparison(object = spata_obj, 
                                 color_by = CD8_exh,
                                 method_gs = "zscore",
                                 smooth = TRUE,
                                 pt_clrsp = "inferno",
                                 smooth_span = 0.2,
                                 pt_size = 2,
                                 display_image = FALSE
)

S1_exh

S1_act <- plotSurfaceComparison(object = spata_obj, 
                                  color_by = CD8_act,
                                  method_gs = "zscore",
                                  smooth = TRUE,
                                  pt_clrsp = "inferno",
                                  smooth_span = 0.2,
                                  pt_size = 2,
                                  display_image = FALSE
)

S1_act


S1_PD <- plotSurfaceComparison(object = spata_obj, 
                                  color_by = PDgene,
                                  method_gs = "zscore",
                                  smooth = TRUE,
                                  pt_clrsp = "inferno",
                                  smooth_span = 0.2,
                                  pt_size = 2,
                                  display_image = FALSE
)

S1_PD

S1_TGF <- plotSurfaceComparison(object = spata_obj, 
                               color_by = TGF,
                               method_gs = "zscore",
                               smooth = TRUE,
                               pt_clrsp = "inferno",
                               smooth_span = 0.2,
                               pt_size = 2,
                               display_image = FALSE
)

S1_TGF

S1_mel <- plotSurfaceComparison(object = spata_obj, 
                                color_by = mel,
                                method_gs = "mean",
                                smooth = TRUE,
                                pt_clrsp = "inferno",
                                smooth_span = 0.2,
                                pt_size = 2,
                                display_image = FALSE
)

S1_mel


S1_GOI <- plotSurfaceComparison(object = spata_obj, 
                                color_by = GOI,
                                method_gs = "mean",
                                smooth = TRUE,
                                pt_clrsp = "inferno",
                                smooth_span = 0.2,
                                pt_size = 2,
                                display_image = FALSE
)

S1_GOI


ggsave(plot = S1_mel, file=file.path("/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATA_object/S4/Gene_Set_Enrichment/zscore/S4_zscore_Melnaoma_markers_Surface_plot2.pdf"),
       width = 6, height = 5)


# plot results
vp <- plotViolinplot(
  object = spata_obj, 
  across = "bayes_space",
  variables = S4mel, 
  clrp = "jama"
)

vp

ggsave(plot = vp, file=file.path("/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATA_object/S4/Gene_Set_Enrichment/zscore/S4_zscore_Mel_markers_violin_plot.pdf"),
       width = 8, height = 6)

# # extract gene set enrichment results
# S1_ssGSEA <- S1_plot[["data"]]
# 
# library(tidyr)
# 
# library(tidyr)
# 
# # Pivot the data frame into wide format
# wide_df <- pivot_wider(S1_ssGSEA, 
#                        id_cols = c(barcodes, x, y, tissue, row, col, imagerow, imagecol, sample, .group),
#                        names_from = variables, 
#                        values_from = values)






#################################################### CYTOSPACE #################################################################

################################# EXTRACT EXPRESSION MATRIX FROM SPATA OBJECT ######################################################

S4mat <- getCountMatrix(S4)
S4cor <- getCoordsDf(S4)
S4cor <- S4cor[, -4:-9]
S4mat[1:5, 1:5]

S4df <- as.data.frame(S4mat)

d <- S4df
names <- rownames(d)
rownames(d) <- NULL
data <- cbind(names,d)
colnames(data)[1] <- "V1"


# Define the file path
file_path <- "S4mat.txt"

# Write data to a tab-delimited text file without row names
write.table(data, file = file_path, sep = "\t", row.names = FALSE, quote = FALSE)

write.table(S4cor, file = "/S4Cor.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# RUN cytospace web tool and get the outputs

# # load library
# -----
library(openxlsx)
library(ggplot2)
library(UCell)
library(gridExtra)
library(reshape2)
library(tidyr)


# -----
# load geneset
# -----
geneset_OF_interest <- read.xlsx("geneset.xlsx")
genesetList <- lapply(1:ncol(geneset_OF_interest), function(i) {geneset_OF_interest[, i]})
names(genesetList) <- colnames(geneset_OF_interest)
genesetList <- lapply(genesetList, function(i) {unique(i[!is.na(i)])})

# -----
# load single cell dataset
# ------- 

melanoma_sc <- read.delim("melanoma_scRNA_GEP.txt", row.names = 1)

# load coordinates

# -----
# get the cells to work
# -----
cells_to_work <- intersect(colnames(melanoma_sc), melanoma_sc_coordinate$OriginalCID)

# -----
# subset single cell dataset
# -----
melanoma_sc <- melanoma_sc[, cells_to_work]


# -----
# deduplicate the coordinates df
# -----
melanoma_sc_coordinate <- dplyr::distinct(melanoma_sc_coordinate, OriginalCID, .keep_all = TRUE)

# -----
# calculate ucell
# -----
melanoma_ucell <- ScoreSignatures_UCell(melanoma_sc, features = genesetList)
melanoma_ucell <- as.data.frame(melanoma_ucell)


# -----
# prepare the df for ggplot
# -----
gg_df <- melanoma_sc_coordinate
# gg_df$LA_tam_ucell <- melanoma_ucell$LA_TAM.Genes_UCell[match(rownames(melanoma_ucell), melanoma_sc_coordinate$OriginalCID)]

# -----
# add ucell scores to the gg dataframe
# -----
melanoma_ucell$cell_name <- rownames(melanoma_ucell)
rownames(melanoma_ucell) <- NULL
gg_df <- merge(gg_df, melanoma_ucell, by.x = "OriginalCID", by.y = "cell_name")
colnames(gg_df) <- gsub("_UCell", "", colnames(gg_df))


# -----
# identify unique spots
# -----
if (FALSE) {
  gg_df$unique_spot <- paste0(gg_df$x, "_", gg_df$y)
  ggplot(gg_df, aes(X, Y, color = unique_spot)) +
    geom_point() +
    theme_bw() + 
    theme(legend.position = "none")  # Remove the legend
}

# -----
# calculate distance
# -----
coordinate_dist <- dist(gg_df[, c("x", "y")], method = "euclidean")
distance_matrix <- as.matrix(coordinate_dist)
min(distance_matrix[distance_matrix != 0])

gg_df$x <- gg_df$x + runif(nrow(gg_df), min = -2, max = 2)
gg_df$y <- gg_df$y + runif(nrow(gg_df), min = -2, max = 2)

# get out dir

outdir <- "/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/Cytoscape_Melanoma/Mel_spatial10x_deconvolution_basedon_Tirosh/results/Cytospace_Output_Ucell/S4"


# # -----
# # plot the enrichment
# # -----
# p <- ggplot(gg_df, aes(x, y, color = CellType)) +
#   geom_point(size = 2) +
#   theme_bw() +
#   scale_colour_brewer(palette = "Paired") +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# p
# 

####
library(ggplot2)
library(dplyr)
library(forcats)
library(patchwork)

# Define the order of CellType for plotting
plot_order <- c("Macrophages", "Melanoma cells", "Fibroblasts","CD8 T cells","CD4 T cells", "B cells",  "Endothelial cells")

# Define custom colors for each CellType
custom_colors <- c("Melanoma cells" = "#F6D776",
                   "Endothelial cells" = "#508D69",
                   "Fibroblasts" = "#FFC5C5",
                   "Macrophages" = "#FF004D",
                   "CD8 T cells" = "#6DB9EF",
                   "CD4 T cells" = "#EEF5FF",
                   "B cells" = "#9ADE7B")

# Reorder gg_df based on plot_order
gg_df <- gg_df %>%
  mutate(CellType = factor(CellType, levels = plot_order))

#plot
p <- ggplot(gg_df, aes(x, y, color = CellType)) +
  geom_point(size = 2) +
  theme_bw() +
  scale_colour_manual(values = custom_colors) +
  guides(alpha = FALSE) +  
  #scale_colour_brewer(palette = "Paired") +
      theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())


print(p)

#save 
ggsave(plot = p, file=file.path(outdir,
                                "cell_type.pdf"),
       width = 8, height = 5)

##
# for (i in names(geneset_OF_interest)) {
#   
#   print(i)
#   p <- ggplot(gg_df, aes(x, y), order = TRUE) +
#     geom_point(aes(color = scale(gg_df[[i]]), order = gg_df[[i]])) +
#     scale_colour_gradient2(low = "lightblue", mid = "lightgrey",  high = "red", na.value = NA) +
#     theme_bw() + ggtitle(i) +
#     labs(color = "UCell score")
#   
#   ggsave(plot = p, file=file.path(outdir,
#                                   paste0(i, ".pdf")),
#          width = 8, height = 5)
#   
# }


#gg_df <- gg_df[, -24]

# try the alternative
library(ggplot2)

# Loop through each gene set
for (i in names(geneset_OF_interest)) {
  print(i)
  
  # Sort gg_df based on the expression of the current feature
  gg_df_sorted <- gg_df[order(gg_df[[i]]), ]
  
  # Create the ggplot object
  p <- ggplot(gg_df_sorted, aes(x, y)) +
    geom_point(aes(color = scale(gg_df_sorted[[i]])), 
               na.rm = TRUE) + # Set na.rm = TRUE to remove NA values
    scale_colour_gradient2(low = "lightblue", mid = "lightgrey",  high = "red", na.value = NA) +
    theme_bw() + 
    ggtitle(i) +
    labs(color = "UCell score")
  
  # Save the plot to a PDF file
  ggsave(plot = p, 
         file = file.path(outdir, paste0(i, ".pdf")),
         width = 8, 
         height = 5)
  
}
  
# -----
# plot the facet ggplot
# -----
new_gg_df <- gg_df
new_gg_df <- melt(new_gg_df, measure.vars = names(geneset_OF_interest))

p <- ggplot(new_gg_df, aes(x, y)) +
  geom_point(aes(color = value)) +
  scale_colour_gradient2(low = "#EEE2DE", high = "#B31312", na.value = NA) +
  theme_bw() + facet_wrap(~variable)
p

# #############
# library(ggplot2)
# library(reshape2)  # for melt function
# 
# # Melt the data
# new_gg_df <- melt(gg_df, measure.vars = names(geneset_OF_interest))
# 
# # Arrange the data frame to have higher values in the front
# new_gg_df <- new_gg_df[order(new_gg_df$value, decreasing = TRUE), ]
# 
# # Create a custom facet grid layout
# grid_layout <- expand.grid(variable = unique(new_gg_df$variable))
# 
# # Define the number of rows and columns for the grid
# n_col <- 3  # number of columns
# n_row <- ceiling(nrow(grid_layout) / n_col)  # number of rows
# 
# # Create plots and store them in a list
# plot_list <- lapply(seq_len(nrow(grid_layout)), function(i) {
#   current_variable <- grid_layout$variable[i]
#   
#   p <- ggplot(subset(new_gg_df, variable == current_variable), aes(x, y)) +
#     geom_point(aes(color = value), size = 0.5) +
#     scale_colour_gradient2(low = "#EEE2DE", high = "#B31312", na.value = NA) +
#     theme_bw() +
#     ggtitle(current_variable) +
#     scale_x_continuous(name = "X Axis", expand = c(0, 0)) + 
#     scale_y_continuous(name = "Y Axis", expand = c(0, 0))
#   
#   return(p)
# })
# 
# # Generate the facet grid plot with different scales
# final_plot <- wrap_plots(plotlist = plot_list, ncol = n_col, scales = "free")
# 
# # Save the final plot to a PDF file
# ggsave(plot = final_plot, 
#        file = file.path(outdir, "facet_grid_plots.pdf"),
#        width = 8, 
#        height = 10)

#######

############ ASSIGN CELL TYPES BASED ON UCELL #################################

# Select the columns of interest
columns_of_interest <- c("IFN_TAMs", "Inflam_TAM", "LA_TAM",
                         "Undifferentiated-Neural.crest-like", "Neural.crest-like", "Neural.crest-like-Transitory",
                         "Transitory", "Transitory-Melanocytic", "Melanocytic",
                         "Exhausion_cell_cyle", "Exhausion_HSP", "Exhausion",
                         "Memory_effector_1", "Early_activated_cells", "Memory_effector_2",
                         "AXL_Melnoma", "MITF_Melanoma")

# Find the column name with the highest value for each row
gg_df$CellTypeFinal <- apply(gg_df[columns_of_interest], 1, function(x) {
  columns_of_interest[which.max(x)]
})

# View the updated gg_df with the new column "CellTypeFinal"
head(gg_df)

# NOW PLOT 

# Define the order of CellType for plotting
plot_order <- c("LA_TAM", "MITF_Melanoma", "Transitory-Melanocytic", "IFN_TAMs", "Inflam_TAM",
                "Undifferentiated-Neural.crest-like", "Neural.crest-like", "Neural.crest-like-Transitory",
                "Transitory", "Melanocytic",
                "Exhausion_cell_cyle", "Exhausion_HSP", "Exhausion",
                "Memory_effector_1", "Early_activated_cells", "Memory_effector_2",
                "AXL_Melnoma")

# Define custom colors for each CellType
custom_colors <- c("MITF_Melanoma" = "#F6D776",
                   "Transitory-Melanocytic" = "#FE7A36",
                   "Undifferentiated-Neural.crest-like" = "#FFF8C9",
                   "Endothelial cells" = "#508D69",
                   "Fibroblasts" = "#FFC5C5",
                   "LA_TAM" = "#FF004D",
                   "CD8 T cells" = "#6DB9EF",
                   "CD4 T cells" = "#EEF5FF",
                   "B cells" = "#9ADE7B")

# Reorder gg_df based on plot_order
gg_df <- gg_df %>%
  mutate(CellTypeFinal = factor(CellTypeFinal, levels = plot_order))

#plot
p <- ggplot(gg_df, aes(x, y, color = CellTypeFinal)) +
  geom_point(size = 2) +
  theme_bw() +
  scale_colour_manual(values = custom_colors) +
  guides(alpha = FALSE) +  
  #scale_colour_brewer(palette = "Paired") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


print(p)

#save 
ggsave(plot = p, file=file.path(outdir,
                                "cell_type.pdf"),
       width = 8, height = 5)


# violin plot 

# plot results
plotViolinplot(
  object = spata_obj, 
  across = "bayes_space",
  variables = TGF, 
  clrp = "jama"
)

#########################################################################################
################################ Infer CNV ##############################################
#########################################################################################



#install iCNV
#BiocManager::install("infercnv")
#laod
library(infercnv)


# rename spata object 

#spata_object <- S1
#spata_object <- S2 
#spata_object <- S3
spata_object <- S4

# see expression matrix

getExpressionMatrixNames(spata_object)

# set data to scaled
setActiveExpressionMatrix(spata_object, "scaled")


# run infer cnv
spata_object <-
  runCnvAnalysis(
    object = spata_object,
    directory_cnv_folder = "7cluster", # example directory
    cnv_prefix = "Chr"
  )

hm <- plotCnvHeatmap(object = spata_object, across = "bayes_space", clrp = "npg")

hm

# surface plot Chr 7 and 10

plotSurface(
  object = spata_object, 
  color_by = "Chr6", 
  pt_clrsp = "Reds",
  pt_size = 3,
  display_image = FALSE, 
  smooth = TRUE, 
  alpha_by = TRUE,
)

plotSurface(
  object = spata_object, 
  color_by = "Chr14", 
  pt_clrsp = "Oslo",
  pt_size = 3,
  display_image = FALSE, 
  smooth = TRUE, 
  alpha_by = TRUE,
)

# compare gene expression on the surface
plot_chr <- plotSurfaceComparison(
  object = spata_object, 
  color_by = c("Chr7", "Chr10"), 
  pt_clrsp = "inferno",
  pt_size = 3,
  display_image = FALSE, 
  smooth = TRUE, 
  alpha_by = TRUE,
) 

plot_chr

##################################################################################################
############################ Segmentation ########################################################
##################################################################################################

#For S1
# create a new segmentation for S1
S4 <- createSpatialSegmentation(S4)

# # interactive plot
plotSurfaceInteractive(object = S4)

# make seg surface plot
S2_seg <- plotSurface(object = S2, color_by = "micro", pt_clrp = "npg", display_image = FALSE)

ggsave(plot = S1_seg, file=file.path("/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATA_object/S1/Segmentation/S1_seg_plot.pdf"),
       width = 6, height = 5)


# vionlin plot across segmentations

# set data to scaled
S1 <- setActiveExpressionMatrix(S1, "scaled")
S1 <- setActiveExpressionMatrix(S1, "denoised")

# get active matrix
getActiveExpressionMatrixName(S1)


# load genes to plot
genes <- c("APOE", "TREM2", "IFNG", "IL10", "IL10RA", "IL10RB", "LDLR", "VLDLR", "LRP1", "LRP5", "LRP8", "APOC1", "C1QB", "CD63", "CD9", "MLANA", "MLANB", "AXL", "MITF")

S1_genes <- c("APOE", "LDLR", "VLDLR", "APOC1", "C1QB", "CD63", "CD9", "MLANA", "MITF")

# names of grouping variables
getGroupingOptions(object = S1)

vp <- plotViolinplot(
  object = S1, 
  variables = S1_genes, 
  across = "micro", 
  ncol = 1, 
  clrp = "npg",
  test_pairwise = "wilcox.test",
  vjust = 1,
  method_gs = "zscore",
  display_facets = TRUE
)

vp









