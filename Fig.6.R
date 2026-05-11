# install necessary packages
BiocManager::install("cytomapper")
BiocManager::install("imcRtools")
devtools::install_github("JinmiaoChenLab/Rphenograph")
BiocManager::install("cytoviewer")
install.packages("Matrix")
BiocManager::install("CATALYST")
install.packages("Matrix")  # Update to the latest version

if (!requireNamespace("BiocParallel", quietly = TRUE)) {
  BiocManager::install("BiocParallel")
}
library(BiocParallel)

# load libraries
library(imcRtools)
library(cytomapper)
library(tidyverse)
library(readxl)
library(dittoSeq)
library(RColorBrewer)
library(Rphenograph)
library(igraph)
library(dittoSeq)
library(viridis)
library(scater)
library(patchwork)
library(cowplot)
library(cytoviewer)
library(pheatmap)

# define path to steinbock
path <- '/Users/...../steinbock-new'
panel_file <- '/Users/..../panel - PENGUIN.csv'

# reading steinbock generated data
IMC <- read_steinbock(path)
IMC

# # load already processed IMC object by Ira
# IMC <- readRDS("/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/NEW_samples_NOV_2024/IMC/MelanomaCohort_IMC/ImagingWorkshop2023-main_Melanoma/data/steinbock-new/R scripts/RDS/CCC_TMA.rds")
# 

# check ounts
counts(IMC)[1:5,1:5]

# check col
head(colData(IMC))

# cehck coords 
head(spatialCoords(IMC))
colPair(IMC, "neighborhood")
head(rowData(IMC))

# Add additional metadata
meta <- read_xlsx('/Users/..../metadata.xlsx')
meta <- meta[, -1]


# # We can set the colnames of the object to generate unique identifiers per cell:
colnames(IMC) <- paste0(IMC$sample_id, "_", IMC$ObjectNumber)
# 
# # add meta data to spatial object 
# # # It is also often the case that sample-CCC_TMAcific metadata are available externally.
# IMC$sample_id = as.character(meta$Patient_ID)
# # It is also often the case that sample-CCC_TMAcific metadata are available externally.
# meta2 <- summary(as.factor(IMC$sample_id)) %>% as.data.frame()
# meta2 <- meta2 %>% dplyr::mutate(ID = rownames(meta2))
# meta2 <- meta2 %>% dplyr::select(ID)

# load R object from Ira
#IMC$indication <- meta$Indication[match(IMC$sample_id, meta$sample_id)]
#CCC_TMA$sample <- meta$Sample_ID[match(CCC_TMA$sample_id, meta$sample_id)]
IMC$ROI <- meta$Tissue_ID[match(IMC$sample_id, meta$sample_id)]
IMC$Response <- meta$Response[match(IMC$sample_id, meta$sample_id)]
IMC$patient_ID <- meta$PIZ[match(IMC$sample_id, meta$sample_id)]

# normalize counts
dittoRidgePlot(IMC, var = "CD8a", group.by = "patient_ID", assay = "counts") +
  ggtitle("CD8a - before transformation")

assay(IMC, "exprs") <- asinh(counts(IMC)/0.5)


dittoRidgePlot(IMC, var = "TREM2", group.by = "patient_ID", split.by= "Response", assay = "exprs") +
  ggtitle("TREM2 - after transformation")

#check APOE and TREM2
dittoRidgePlot(IMC, var = "TREM2", group.by = "Response", assay = "exprs") +
  ggtitle("TREM2")

#save IMC object
saveRDS(IMC, '/Users/..../IMC.rds')

# load IMC object

IMC <- readRDS('/Users/...../IMC.rds')

################################ Phenoclustering #########################################################################
# Example of using this palette in a plot

color_vectors <- list()

my_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
               "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
               "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", "#393b79",
               "#637939", "#8c6d31", "#843c39", "#7b4173", "#3182bd", "#6baed6", "#9e9ac8",
               "#fd8d3c", "#e6550d", "#31a354", "#756bb1", "#636363", "#e41a1c", "#377eb8",
               "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999",
               "#fb8072", "#80b1d3", "#b3de69")

ROI <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
         "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
         "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", "#393b79",
         "#637939", "#8c6d31", "#843c39", "#7b4173", "#3182bd", "#6baed6", "#9e9ac8",
         "#fd8d3c", "#e6550d", "#31a354", "#756bb1", "#636363", "#e41a1c", "#377eb8",
         "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999",
         "#fb8072", "#80b1d3", "#b3de69")

sample <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
            "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
            "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", "#393b79",
            "#637939", "#8c6d31", "#843c39", "#7b4173", "#3182bd", "#6baed6", "#9e9ac8",
            "#fd8d3c", "#e6550d", "#31a354", "#756bb1", "#636363", "#e41a1c", "#377eb8",
            "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999",
            "#fb8072", "#80b1d3", "#b3de69")

Response <- c("#d6604d", "#4393c3")

Patient_ID <- c("#ff7f0e", "#2ca02c", "#d62728", "#ffff33", "#f781bf", "#3182bd")

color_vectors$ROI <- ROI
color_vectors$sample <- sample
color_vectors$Response <- Response
color_vectors$Patient_ID <- Patient_ID
#color_vectors$immunetype <- immuntype
metadata(IMC)$color_vectors <- color_vectors

#for heat maps
rowData(IMC)$use_channel_heat <- !grepl("DNA1|DNA2|ki67|Hoechst|HH3|Collagen1|SMA|HLA-DR|BCatenin|cCasp3|CytC", rownames(IMC))
#for umap
rowData(IMC)$use_channel <- !grepl("DNA1|DNA2|ki67|Hoechst|HH3|Collagen1|SMA|HLA-DR|BCatenin|cCasp3|CytC", rownames(IMC))

cur_cells <- sample(seq_len(ncol(IMC)), 20000)

# get expression matrix
mat <- t(assay(IMC, "exprs")[rowData(IMC)$use_channel,])

# #scale rowwise
# scaled_mat <- t(scale(t(mat)))
# any(is.na(scaled_mat))   # Check for NA
# any(is.nan(scaled_mat))  # Check for NaN
# any(is.infinite(scaled_mat))  # Check for Inf
# which(is.na(scaled_mat), arr.ind = TRUE)   # Location of NAs
# which(is.nan(scaled_mat), arr.ind = TRUE)  # Location of NaNs
# scaled_mat[is.na(scaled_mat)] <- 0
# scaled_mat[is.nan(scaled_mat)] <- 0

# run phenograph
set.seed(230619)
out <- Rphenograph(mat, k = 145)

# Convert the clustering results into a vector
clusters <- factor(membership(out[[2]]))

# add a new column in IMC object 
IMC$clusters <- clusters

# Define a custom color palette
custom_colors <- colorRampPalette(c("blue", "white", "red"))(50)

# Set breaks for the heatmap to range from -2 to +2
breaks_seq <- seq(-2, 2, length.out = 50)


hm <- dittoHeatmap(IMC[,cur_cells], 
                   genes = rownames(IMC) [rowData(IMC)$use_channel],
                   assay = "exprs", scale = "row",
                   heatmap.colors = custom_colors,  # Custom color palette
                   breaks = breaks_seq,  # Breaks from -2 to +2
                   annot.by = c("clusters", "patient_ID", "Response"),
                   annot.colors = c(dittoColors(1)[1:length(unique(IMC$clusters))],
                                    metadata(IMC)$color_vectors$sample,
                                    metadata(IMC)$color_vectors$Response,
                                    metadata(IMC)$color_vectors$Patient_ID))


hm

dittoDimPlot(IMC, var = "pg_clusters", 
             reduction.use = "UMAP", size = 0.2,
             do.label = TRUE) +
  ggtitle("Phenograph clusters on UMAP")

######################### seubsetting the IMC spatial object #################################################
unique(IMC$Response)

# subset based on response column 
IMC_non_responder <- IMC[, IMC$Response == "Non-Responder"] 
IMC_responder <- IMC[, IMC$Response == "Responder"]

#check metadata of the subseted object 
rownames(IMC_responder)

dittoRidgePlot(IMC_responder, var = "CD8a", group.by = "patient_ID", assay = "exprs") +
  ggtitle("CD8a Expression in Responders")

dittoRidgePlot(IMC_non_responder, var = "CD8a", group.by = "Response", assay = "exprs") +
  ggtitle("CD8a Expression in Non-Responders")


# non-responder 

#for heat maps
rowData(IMC_non_responder)$use_channel_heat <- !grepl("DNA1|DNA2|ki-67|Hoechst|HH3|Collagen|SMA|HLA-DR|BCatenin|FoxP3|cCasp3|CytC|PanKer|CD11b|Lag3", rownames(IMC_non_responder))
#for umap
rowData(IMC_non_responder)$use_channel <- !grepl("DNA1|DNA2|ki-67|Hoechst|HH3|Collagen|SMA|HLA-DR|BCatenin|FoxP3|cCasp3|CytC|PanKer|CD11b|Lag3", rownames(IMC_non_responder))

set.seed(1234)
cur_cells <- sample(seq_len(ncol(IMC_non_responder)), 10000)

# get expression matrix
mat_non_responder <- t(assay(IMC_non_responder, "exprs")[rowData(IMC_non_responder)$use_channel,])

set.seed(230619)
out <- Rphenograph(mat_non_responder, k = 195)

clusters <- factor(membership(out[[2]]))

IMC_non_responder$pg_clusters <- clusters

# Define a custom color palette
custom_colors <- colorRampPalette(c("blue", "white", "red"))(50)

# Set breaks for the heatmap to range from -2 to +2
breaks_seq <- seq(-2, 2, length.out = 50)

hm <- dittoHeatmap(IMC_non_responder[,cur_cells], 
                   genes = rownames(IMC_non_responder) [rowData(IMC_non_responder)$use_channel],
                   assay = "exprs", scale = "row",
                   annotation_colors = viridis(100),
                   heatmap.colors = custom_colors,  # Custom color palette
                   breaks = breaks_seq,  # Breaks from -2 to +2
                   annot.by = c("pg_clusters", "patient_ID", "Response", "ROI"),
                   annot.colors = c(dittoColors(1)[1:length(unique(IMC_non_responder$pg_clusters))],
                                    metadata(IMC_non_responder)$color_vectors$sample,
                                    metadata(IMC_non_responder)$color_vectors$Response,
                                    metadata(IMC_non_responder)$color_vectors$Patient_ID,
                                    metadata(IMC_non_responder)$color_vectors$ROI))


hm

library(dplyr)
cluster_celltype <- recode(IMC_non_responder$pg_clusters,
                           "1" = "VIM_MES",
                           "2" = "CD8",
                           "3" = "Macrophage_ALA",
                           "4" = "VIM_MES",
                           "5" = "Tumor",
                           "6" = "Macrophage",
                           "7" = "Tumor_MES",
                           "8" = "Memory_effector",
                           "9" = "UK",
                           "10" = "Bcell",
                           "11" = "Vascular_Macrophage_IFN",
                           "12" = "Tumor",
                           "13" = "HSC",
                           "14" = "Memory_CD45RO",
                           '15' = "HSC",
                           "16" = "Proliferating_cell",
                           "17" = "UK",
                           "18" = "Macrophage_CD38",
                           "19" = "UK",
                           "20" = "Tumor_stemcell",
                           "21" = "Proliferating_cell",
                           "22" = "UK",
                           "23" = "Endothelial_cell",
                           "24" = "UK",
                           "25" = "Tumor",
                           "26" = "Macrophage",
                           "27" = "UK",
                           "28" = "Tumor",
                           "29"= "Proliferating_tumor")

IMC_non_responder$cluster_celltype <- cluster_celltype

# now create a broader cell type category

celltype_braod <- recode(IMC_non_responder$cluster_celltype,
                         "VIM_MES" = "MES",
                         "CD8" = "CD8",
                         "Macrophage_ALA" = "Macrophage",
                         "VIM_MES" = "MES",
                         "Tumor" = "Tumor",
                         "Macrophage" = "Macrophage",
                         "Tumor_MES" = "Tumor",
                         "Memory_effector" = "Memory",
                         "UK" = "UK",
                         "Bcell" = "Bcell",
                         "Vascular_Macrophage_IFN" = "Macrophage",
                         "Tumor" = "Tumor",
                         "HSC" = "HSC",
                         "Memory_CD45RO" = "Memory",
                         'HSC' = "HSC",
                         "Proliferating_cell" = "Proliferating_cell",
                         "UK" = "UK",
                         "Macrophage_CD38" = "Macrophage",
                         "UK" = "UK",
                         "Tumor_stemcell" = "Tumor",
                         "Proliferating_cell" = "Proliferating_cell",
                         "UK" = "UK",
                         "Endothelial_cell" = "Endothelial_cell",
                         "UK" = "UK",
                         "Tumor" = "Tumor",
                         "Macrophage" = "Macrophage",
                         "UK" = "UK",
                         "Tumor" = "Tumor",
                         "Proliferating_tumor"= "Tumor")

IMC_non_responder$celltype_braod <- celltype_braod

# save non-responder IMC object 
saveRDS(IMC_non_responder, '/Users/...../non_Responder_IMC.rds')

########################################## Images Processing and Visualization ######################################

# define image path
img_path <- '/Users/...../steinbock-new/img'

#load images
images <- loadImages(img_path)

# set mask
masks <- loadImages('/Users/...../mesmer_masks', as.is = TRUE)

#set chanel names
channelNames(images) <- rownames(IMC)
images

#Here, we will store the matched sample_id, patient_id and indication information within the elementMetadata slot of the multi-channel images and segmentation masks objects. It is crucial that the order of the images in both CytoImageList objects is the same.
all.equal(names(images), names(masks))

# Extract patient id from image name
patient_id <- names(images)

# Retrieve cancer type per patient from metadata file
Response <- meta$Response[match(patient_id, meta$sample_id)] 

# Store patient and image level information in elementMetadata
mcols(images) <- mcols(masks) <- DataFrame(sample_id = names(images),
                                           patient_id = patient_id,
                                           Response = Response)


#Cytomapper
cytomapper_IMC <- measureObjects(masks, image = images, img_id = "sample_id")

cytomapper_IMC

#Image visualization
library(SpatialExperiment)
library(cytomapper)

IMC <- readRDS('/Users/../IMC.rds')
images <- readRDS('/Users/../images.rds')
masks <- readRDS('/Users/../masks.rds')

# Sample images
set.seed(220517)
# Predefined IDs you want to include

# non_res_id
non_res_id <- c(
  "28_1_1", "28_2_2", "28_3_3",
  "20_1_1", "20_2_2", "20_3_3", "20_4_4",
  "09_1_1", "09_2_2", "09_3_3", "09_4_4"
)

# Remaining IDs (excluding the specific ones)
res_id <- setdiff(unique(IMC$sample_id), non_res_id)

# For Non Responder Patients
#cur_id <- sample(unique(IMC$sample_id), 9)
cur_images <- images[names(images) %in% non_res_id] # non_responder patients
cur_masks <- masks[names(masks) %in% non_res_id] # non_responder patients

# For Responder Patients
cur_images <- images[names(images) %in% res_id] # responder patients
cur_masks <- masks[names(masks) %in% res_id] # responder patients

#Pixel visualization APOE+ TERM2
plotPixels(cur_images, 
           colour_by = c("S100","CD45", "CD31"),
           bcg = list(S100 = c(0, 2, 1),
                      CD45 = c(0, 5, 1),
                      CD31 = c(0, 5, 1)))


#Pixel visualization APOE+ TERM2
plotPixels(cur_images, 
           colour_by = c("TREM2","CD8a", "GranzymeB"),
           bcg = list(TREM2 = c(0, 10, 1),
                      CD8a = c(0, 5, 1),
                      GranzymeB = c(0, 10, 1)))


# ONLY TREM2
#Pixel visualization
plotPixels(cur_images, 
           colour_by = "HLADR",
           bcg = list(HLADR = c(0, 5, 1)),
           colour = list(HLADR = c("black", "cyan2")))


# visualize cells
plotCells(cur_masks,
          object = IMC, 
          cell_id = "ObjectNumber", 
          img_id = "sample_id",
          colour_by = c("CD8a", "TREM2", "CD163"),
          exprs_values = "exprs",
          colour = list(CD8a = c("black", "burlywood1"),
                        TREM2 = c("black", "cyan2"),
                        CD163 = c("black", "firebrick1")))


# define cell colors
# Get unique levels
celltype_levels <- unique(IMC_non_responder$celltype_braod)

# Define colors for each level
color_celltype_braod <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
  "#8c564b", "#e377c2", "#ffff33", "#999999", "#003049", "#8dd3c7"
)

# Verify lengths match
if (length(color_celltype_braod) < length(celltype_levels)) {
  stop("Not enough colors for celltype levels.")
}

# Assign colors dynamically
names(color_celltype_braod) <- celltype_levels
colour_list <- list(celltype_braod = color_celltype_braod)

plotPixels(
  image = cur_images,
  mask = cur_masks,
  object = IMC_non_responder, 
  cell_id = "ObjectNumber", 
  img_id = "sample_id",
  colour_by = c("CD8a", "TREM2", "CD163"),
  outline_by = "celltype_braod",
  bcg = list(CD8a = c(0, 5, 1), TREM2 = c(0, 5, 1), CD163 = c(0, 5, 1)),
  colour = colour_list,  # Fixed color mapping
  thick = TRUE
)


# CYTOVIEWER
library(cytoviewer)

app <- cytoviewer(image = images,
                  mask = masks,
                  object = IMC,
                  cell_id = "ObjectNumber",
                  img_id = "sample_id")

if (interactive()) {
  shiny::runApp(app)
}

######################################### REMOVING UNINFORMATIVE/EMPTY/ ROIS #########################
unique(IMC$ROI)
unique(IMC$sample_id)

# subset based on response column 
IMC_filtered <- IMC[, !(IMC$sample_id %in% IDs_with_low_quality)]

#check metadata of the subseted object 
rownames(IMC_filtered)

################################ Phenoclustering on Filtered object #########################################################################
# Example of using this palette in a plot

color_vectors <- list()

my_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
               "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
               "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", "#393b79",
               "#637939", "#8c6d31", "#843c39", "#7b4173", "#3182bd", "#6baed6", "#9e9ac8",
               "#fd8d3c", "#e6550d", "#31a354", "#756bb1", "#636363", "#e41a1c", "#377eb8",
               "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999",
               "#fb8072", "#80b1d3", "#b3de69")

ROI <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
         "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
         "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", "#393b79",
         "#637939", "#8c6d31", "#843c39", "#7b4173", "#3182bd", "#6baed6", "#9e9ac8",
         "#fd8d3c", "#e6550d", "#31a354", "#756bb1", "#636363", "#e41a1c", "#377eb8",
         "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999",
         "#fb8072", "#80b1d3", "#b3de69")

sample <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
            "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
            "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", "#393b79",
            "#637939", "#8c6d31", "#843c39", "#7b4173", "#3182bd", "#6baed6", "#9e9ac8",
            "#fd8d3c", "#e6550d", "#31a354", "#756bb1", "#636363", "#e41a1c", "#377eb8",
            "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999",
            "#fb8072", "#80b1d3", "#b3de69")

Response <- c("#d6604d", "#4393c3")

Patient_ID <- c("#ff7f0e", "#2ca02c", "#d62728", "#ffff33", "#f781bf", "#3182bd")

celltype_broad <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2")

cluster_celltype <- c("#ff7f0e", "#2ca02c", "#d62728", "#ffff33", "#f781bf", "#3182bd",
                      "#843c39", "#fb8072", "#b3de69", "#393b79", "#8c6d31")

color_vectors$ROI <- ROI
color_vectors$sample <- sample
color_vectors$Response <- Response
color_vectors$Patient_ID <- Patient_ID
color_vectors$celltype_braod <- celltype_broad
color_vectors$cluster_celltype <- cluster_celltype
metadata(IMC_filtered)$color_vectors <- color_vectors

#for heat maps
rowData(IMC_filtered)$use_channel_heat <- !grepl("DNA1|DNA2|ki67|Hoechst|HH3|Collagen1|SMA|HLA-DR|BCatenin|cCasp3|CytC", rownames(IMC_filtered))
#for umap
rowData(IMC_filtered)$use_channel <- !grepl("DNA1|DNA2|ki67|Hoechst|HH3|Collagen1|SMA|HLA-DR|BCatenin|cCasp3|CytC", rownames(IMC_filtered))

cur_cells <- sample(seq_len(ncol(IMC_filtered)), 20000)

# get expression matrix
mat <- t(assay(IMC_filtered, "exprs")[rowData(IMC_filtered)$use_channel,])

# #scale rowwise
# scaled_mat <- t(scale(t(mat)))
# any(is.na(scaled_mat))   # Check for NA
# any(is.nan(scaled_mat))  # Check for NaN
# any(is.infinite(scaled_mat))  # Check for Inf
# which(is.na(scaled_mat), arr.ind = TRUE)   # Location of NAs
# which(is.nan(scaled_mat), arr.ind = TRUE)  # Location of NaNs
# scaled_mat[is.na(scaled_mat)] <- 0
# scaled_mat[is.nan(scaled_mat)] <- 0

# run phenograph
set.seed(230619)
out <- Rphenograph(mat, k = 145)

# Convert the clustering results into a vector
clusters <- factor(membership(out[[2]]))

# add a new column in IMC object 
IMC_filtered$clusters <- clusters

# Define a custom color palette
custom_colors <- colorRampPalette(c("blue", "white", "red"))(50)

# Set breaks for the heatmap to range from -2 to +2
breaks_seq <- seq(-1, 1, length.out = 50)


hm <- dittoHeatmap(IMC_filtered[,cur_cells], 
                   genes = rownames(IMC_filtered) [rowData(IMC_filtered)$use_channel],
                   assay = "exprs", scale = "row",
                   heatmap.colors = custom_colors,  # Custom color palette
                   breaks = breaks_seq,  # Breaks from -2 to +2
                   annot.by = c("clusters", "patient_ID", "Response"),
                   annot.colors = c(dittoColors(1)[1:length(unique(IMC_filtered$clusters))],
                                    metadata(IMC_filtered)$color_vectors$sample,
                                    metadata(IMC_filtered)$color_vectors$Response,
                                    metadata(IMC_filtered)$color_vectors$Patient_ID))


hm

library(dplyr)
cluster_celltype <- recode(IMC_filtered$clusters,
                           "1" = "Macrophage",
                           "2" = "Tcell_CD8",
                           "3" = "Macrophage_ALA",
                           "4" = "UK",
                           "5" = "HSC",
                           "6" = "Endothelial_cell",
                           "7" = "Endothelial_cell",
                           "8" = "Macrophage_IFN",
                           "9" = "Macrophage",
                           "10" = "Tumor",
                           "11" = "Tumor",
                           "12" = "Tumor",
                           "13" = "UK",
                           "14" = "Tumor",
                           '15' = "Tcell_Memory_Effector",
                           "16" = "Macrophage_IFN",
                           "17" = "UK",
                           "18" = "UK",
                           "19" = "Tumor",
                           "20" = "Tumor",
                           "21" = "Tumor",
                           "22" = "UK",
                           "23" = "UK",
                           "24" = "Tcell_Effector",
                           "25" = "UK",
                           "26" = "UK",
                           "27" = "Bcell",
                           "28" = "UK",
                           "29" = "HSC",
                           "30" = "Tumor",
                           "31" = "UK",
                           "32" = "Tumor",
                           "33" = "Endothelial_cell")

IMC_filtered$cluster_celltype <- cluster_celltype

# now create a broader cell type category
celltype_broad <- recode(IMC_filtered$cluster_celltype,
                         "Macrophage" = "Macrophage",
                         "Tcell_CD8" = "Tcell",
                         "Macrophage_ALA" = "Macrophage",
                         "UK" = "Unknown",
                         "HSC" = "HSC" ,
                         "Endothelial_cell" = "Endothelial",
                         "Endothelial_cell" = "Endothelial",
                         "Macrophage_IFN" = "Macrophage",
                         "Macrophage" = "Macrophage",
                         "Tumor" = "Tumor",
                         "Tumor" = "Tumor",
                         "Tumor" = "Tumor",
                         "UK" = "Unknown",
                         "Tumor" = "Tumor",
                         "Tcell_Memory_Effector" = "Tcell",
                         'Macrophage_IFN' = "Macrophage",
                         "UK" = "Unknown",
                         "UK" = "Unknown",
                         "Tumor" = "Tumor",
                         "Tumor" = "Tumor",
                         "Tumor" = "Tumor",
                         "UK" = "Unknown",
                         "UK" = "UK",
                         "Bcell" = "Bcell",
                         "Tcell_Effector" = "Tcell",
                         "UK" = "Unknown",
                         "HSC" = "HSC",
                         "Tumor" = "Tumor",
                         "UK" = "Unknown",
                         "Tumor" = "Tumor",
                         "Endothelial_cell"= "Endothelial")

IMC_filtered$celltype_broad <- celltype_broad

# remove unknown cells
# Check the levels of celltype_broad
unique(IMC_filtered$celltype_broad)

# Subset the spatial object to exclude "Unknown"
IMC_filtered <- IMC_filtered[, IMC_filtered$celltype_broad != "Unknown"]

# remove wrong column 
IMC_filtered$celltype_braod <- NULL

# Verify the subset worked
unique(IMC_filtered$celltype_broad)

# check rownames
rownames(IMC_filtered)

# heatmap cluster vs. cell type
tab1 <- table(IMC_filtered$celltype_broad, 
              paste("Rphenograph", IMC_filtered$clusters))

tab2 <- table(IMC_filtered$cluster_celltype, 
              paste("Rphenograph", IMC_filtered$clusters))

pheatmap(log10(tab1 + 10), color = viridis(100))
pheatmap(log10(tab2 + 10), color = viridis(100))

set.seed(220818)
cur_cells <- sample(seq_len(ncol(IMC_filtered)), 20000)

# Select genes based on the updated marker_class
genes_detailed_celltype <- c("CD163", "CD204", "CD68", "APOE", "TREM2", "CCL2", "CCL4", "S100", "CD8a", "CD45RO", "CD3", "CD34", "CD31", "CD20", "GranzymeB", "PDL1", "PD1", "Tim3", "Lag3", "TIGIT", "Tbet")

gene_braod_celltype <- c("CD163", "CD204", "CD68", "S100", "CD8a", "CD3", "CD34", "CD31", "CD20")

GOI <- c("APOE", "TREM2", "CCL2", "CCL4") 


# Convert Response to a factor
colData(IMC_filtered)$Response <- factor(colData(IMC_filtered)$Response, 
                                         levels = c("Responder", "Non-Responder"))

# Define a custom color palette
custom_colors <- colorRampPalette(c("blue", "white", "red"))(50)

# Set breaks for the heatmap to range from -2 to +2
breaks_seq <- seq(-1, 1, length.out = 50)

####################################### Detailed Clusters ########################################################
# plot heatmap all genes
hm <- dittoHeatmap(IMC_filtered[,cur_cells], 
                   genes = rownames(IMC_filtered) [rowData(IMC_filtered)$use_channel],
                   assay = "exprs", scale = "row",
                   heatmap.colors = custom_colors,  # Custom color palette
                   breaks = breaks_seq,  # Breaks from -2 to +2
                   annot.by = c("cluster_celltype", "patient_ID", "Response"),
                   annot.colors = c(dittoColors(1)[1:length(unique(IMC_filtered$clusters))],
                                    metadata(IMC_filtered)$color_vectors$cluster_celltype,
                                    metadata(IMC_filtered)$color_vectors$Response,
                                    metadata(IMC_filtered)$color_vectors$Patient_ID))


hm

# plot heatmap only marker genes
library(scuttle)

# Use aggregateAcrossCells with proper metadata propagation
celltype_mean <- aggregateAcrossCells(
  as(IMC_filtered, "SingleCellExperiment"), 
  ids = IMC_filtered$cluster_celltype,  # Group by cell type
  statistics = "mean",
  use.assay.type = "exprs", 
  subset.row = genes_detailed_celltype
)

# Create the heatmap
dittoHeatmap(
  celltype_mean,
  assay = "exprs", 
  cluster_cols = TRUE, 
  scale = "row",
  heatmap.colors = custom_colors,  # Custom color palette
  annot.by = "cluster_celltype"  # Ensure these exist in colDat# Corrected annotation colors
)

####################################### Broad Clusters ########################################################
# plot heatmap all genes
hm <- dittoHeatmap(IMC_filtered[,cur_cells], 
                   genes = rownames(IMC_filtered) [rowData(IMC_filtered)$use_channel],
                   assay = "exprs", scale = "row",
                   heatmap.colors = custom_colors,  # Custom color palette
                   breaks = breaks_seq,  # Breaks from -2 to +2
                   annot.by = c("celltype_broad", "patient_ID", "Response"),
                   annot.colors = c(dittoColors(1)[1:length(unique(IMC_filtered$clusters))],
                                    metadata(IMC_filtered)$color_vectors$cluster_celltype,
                                    metadata(IMC_filtered)$color_vectors$Response,
                                    metadata(IMC_filtered)$color_vectors$Patient_ID))


hm

# plot heatmap only marker genes
library(scuttle)

# Use aggregateAcrossCells with proper metadata propagation
celltype_mean <- aggregateAcrossCells(
  as(IMC_filtered, "SingleCellExperiment"), 
  ids = IMC_filtered$celltype_broad,  # Group by cell type
  statistics = "mean",
  use.assay.type = "exprs", 
  subset.row = genes_detailed_celltype
)

# Create the heatmap
dittoHeatmap(
  celltype_mean,
  assay = "exprs", 
  cluster_cols = TRUE, 
  scale = "row",
  heatmap.colors = custom_colors,  # Custom color palette
  annot.by = "celltype_broad"  # Ensure these exist in colDat# Corrected annotation colors
)

############################# subsetting by cell type and response ########################################

######################### seubsetting the IMC_filtered spatial object #################################################
# load IMC object if necessary
IMC_filtered <- readRDS('/Users/..../IMC_Filtered_UK_removed.rds')

unique(IMC_filtered$Response)

# subset based on response column 
IMC_non_responder <- IMC_filtered[, IMC_filtered$Response == "Non-Responder"] 
IMC_responder <- IMC_filtered[, IMC_filtered$Response == "Responder"]

#check metadata of the subseted object 
rownames(IMC_responder)

# Define a custom color palette
custom_colors <- colorRampPalette(c("blue", "white", "red"))(50)

# Set breaks for the heatmap to range from -2 to +2
breaks_seq <- seq(-1, 1, length.out = 50)

##################################### non-responder analysis ###########################################################

set.seed(220818)
cur_cells <- sample(seq_len(ncol(IMC_non_responder)), 20000)

# plot heatmap all genes
hm <- dittoHeatmap(IMC_non_responder[,cur_cells], 
                   genes = rownames(IMC_non_responder) [rowData(IMC_non_responder)$use_channel],
                   assay = "exprs", scale = "row",
                   heatmap.colors = custom_colors,  # Custom color palette
                   breaks = breaks_seq,  # Breaks from -2 to +2
                   annot.by = c("cluster_celltype", "patient_ID", "Response"),
                   annot.colors = c(dittoColors(1)[1:length(unique(IMC_non_responder$clusters))],
                                    metadata(IMC_non_responder)$color_vectors$cluster_celltype,
                                    metadata(IMC_non_responder)$color_vectors$Response,
                                    metadata(IMC_non_responder)$color_vectors$Patient_ID))


hm

# plot heatmap only marker genes
library(scuttle)

# Use aggregateAcrossCells with proper metadata propagation
celltype_mean <- aggregateAcrossCells(
  as(IMC_non_responder, "SingleCellExperiment"), 
  ids = IMC_non_responder$cluster_celltype,  # Group by cell type
  statistics = "mean",
  use.assay.type = "exprs", 
  subset.row = genes_detailed_celltype
)

# Create the heatmap
dittoHeatmap(
  celltype_mean,
  assay = "exprs", 
  cluster_cols = TRUE, 
  scale = "row",
  scaled.to.max = FALSE,
  heatmap.colors = custom_colors,  # Custom color palette
  annot.by = "cluster_celltype"  # Ensure these exist in colDat# Corrected annotation colors
)

# stacked bar plot non-responder
# by sample_id - percentage
dittoBarPlot(IMC_non_responder, 
             var = "cluster_celltype", 
             group.by = "sample_id") +
  scale_fill_manual(values = metadata(IMC_non_responder)$color_vectors$cluster_celltype)


##################################### responder analysis ###########################################################

set.seed(220818)
cur_cells <- sample(seq_len(ncol(IMC_responder)), 20000)

# plot heatmap all genes
hm <- dittoHeatmap(IMC_responder[,cur_cells], 
                   genes = rownames(IMC_responder) [rowData(IMC_responder)$use_channel],
                   assay = "exprs", scale = "row",
                   heatmap.colors = custom_colors,  # Custom color palette
                   breaks = breaks_seq,  # Breaks from -2 to +2
                   annot.by = c("cluster_celltype", "patient_ID", "Response"),
                   annot.colors = c(dittoColors(1)[1:length(unique(IMC_responder$clusters))],
                                    metadata(IMC_responder)$color_vectors$cluster_celltype,
                                    metadata(IMC_responder)$color_vectors$Response,
                                    metadata(IMC_responder)$color_vectors$Patient_ID))


hm

# plot heatmap only marker genes
library(scuttle)

# Use aggregateAcrossCells with proper metadata propagation
celltype_mean <- aggregateAcrossCells(
  as(IMC_responder, "SingleCellExperiment"), 
  ids = IMC_responder$cluster_celltype,  # Group by cell type
  statistics = "mean",
  use.assay.type = "exprs", 
  subset.row = genes_detailed_celltype
)

# Create the heatmap
dittoHeatmap(
  celltype_mean,
  assay = "exprs", 
  cluster_cols = TRUE, 
  scale = "row",
  scaled.to.max = FALSE,
  heatmap.colors = custom_colors,  # Custom color palette
  annot.by = "cluster_celltype"  # Ensure these exist in colDat# Corrected annotation colors
)

# Violin Plot - plotExpression
plotExpression(IMC_responder[,cur_cells], 
               features = genes_detailed_celltype,
               x = "cluster_celltype", 
               exprs_values = "exprs", 
               colour_by = "cluster_celltype") +
  theme(axis.text.x =  element_text(angle = 90))+
  scale_color_manual(values = metadata(IMC_responder)$color_vectors$cluster_celltype)

###################################### cell type proportion in responder and non-responder patients #################################################
# load IMC object if necessary
IMC_filtered <- readRDS('/Users/.../IMC_Filtered.rds')

# subset based on response column 
IMC_non_responder <- IMC_filtered[, IMC_filtered$Response == "Non-Responder"] 
IMC_responder <- IMC_filtered[, IMC_filtered$Response == "Responder"]

# by sample_id - percentage nonresponder
non_res <- dittoBarPlot(IMC_non_responder, 
                        var = "cluster_celltype", 
                        group.by = "sample_id") +
  scale_fill_manual(values = metadata(IMC_non_responder)$color_vectors$cluster_celltype)

# by sample_id - percentage responder
res <- dittoBarPlot(IMC_responder, 
                    var = "cluster_celltype", 
                    group.by = "sample_id") +
  scale_fill_manual(values = metadata(IMC_responder)$color_vectors$cluster_celltype)

########################################## agreegate bar plots ##########################################
library(dplyr)
library(ggplot2)

# Extract relevant data for aggregation
non_responder_data <- as.data.frame(colData(IMC_non_responder))  # Extract column data
responder_data <- as.data.frame(colData(IMC_responder))          # Extract column data

# Add a "Response" column to both datasets
non_responder_data <- non_responder_data %>%
  mutate(Response = "Non-Responder")
responder_data <- responder_data %>%
  mutate(Response = "Responder")

# Combine the two datasets
combined_data <- bind_rows(non_responder_data, responder_data)

# Calculate cell type proportions within each Response group
stacked_data <- combined_data %>%
  group_by(Response, cluster_celltype) %>%
  summarize(cell_count = n(), .groups = "drop") %>%
  group_by(Response) %>%
  mutate(average_percentage = cell_count / sum(cell_count) * 100)

# Create the stacked bar plot
stacked_bar_plot <- ggplot(stacked_data, aes(x = Response, y = average_percentage, fill = cluster_celltype)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = metadata(IMC_filtered)$color_vectors$cluster_celltype) +
  theme_minimal() +
  labs(title = "Cell Type Proportions in Responders vs. Non-Responders",
       x = "Response Category", y = "Average Percentage", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display the plot
print(stacked_bar_plot)

######################################### t-test and zitter plot #########################################
library(dplyr)
library(ggplot2)

# Extract relevant data
data <- as.data.frame(colData(IMC_filtered))

# Calculate percentages for each sample and cell type
sample_data <- data %>%
  group_by(sample_id, Response, cluster_celltype) %>%
  summarize(cell_count = n(), .groups = "drop") %>%
  group_by(sample_id, Response) %>%
  mutate(percentage = cell_count / sum(cell_count) * 100)

# Filter data for Macrophage_ALA and Macrophage_IFN
filtered_data <- sample_data %>%
  filter(cluster_celltype %in% c("Macrophage_ALA", "Macrophage_IFN"))

# Perform t-tests for Macrophage_ALA and Macrophage_IFN
t_test_ala <- t.test(percentage ~ Response, data = filtered_data %>% filter(cluster_celltype == "Macrophage_ALA"))
t_test_ifn <- t.test(percentage ~ Response, data = filtered_data %>% filter(cluster_celltype == "Macrophage_IFN"))

# Print t-test results
cat("T-Test for Macrophage_ALA:\n")
print(t_test_ala)
cat("\nT-Test for Macrophage_IFN:\n")
print(t_test_ifn)

# Create jittered box plot
jitter_box_plot <- ggplot(filtered_data, aes(x = Response, y = percentage, fill = cluster_celltype)) +
  geom_boxplot(alpha = 0.5) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, aes(color = cluster_celltype)) +
  facet_wrap(~cluster_celltype) +
  scale_fill_manual(values = c("Macrophage_ALA" = "blue", "Macrophage_IFN" = "red")) +
  scale_color_manual(values = c("Macrophage_ALA" = "blue", "Macrophage_IFN" = "red")) +
  theme_minimal() +
  labs(title = "Macrophage_ALA and Macrophage_IFN Percentages by Response",
       x = "Response Category", y = "Percentage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  annotate("text", x = 1.5, y = max(filtered_data$percentage) + 5,
           label = sprintf("ALA p=%.3f\nIFN p=%.3f", t_test_ala$p.value, t_test_ifn$p.value),
           color = "black", size = 5, hjust = 0.5)

# Display the plot
print(jitter_box_plot)

########################################## Images Processing and Visualization ######################################
#load meta data
# Add additional metadata
meta <- read_xlsx('/Users/../metadata.xlsx')
meta <- meta[, -1]

# load IMC object with with all ROIS
IMC <- readRDS("/Users/../IMC.rds")

# define image path
img_path <- '/Users/../steinbock-new/img'

#load images
images <- loadImages(img_path)

# set mask
masks <- loadImages('/Users../mesmer_masks', as.is = TRUE)

#set chanel names
channelNames(images) <- rownames(IMC)
images

#Here, we will store the matched sample_id, patient_id and indication information within the elementMetadata slot of the multi-channel images and segmentation masks objects. It is crucial that the order of the images in both CytoImageList objects is the same.
all.equal(names(images), names(masks))

# Extract patient id from image name
patient_id <- names(images)

# Retrieve cancer type per patient from metadata file
Response <- meta$Response[match(patient_id, meta$sample_id)] 

# Store patient and image level information in elementMetadata
mcols(images) <- mcols(masks) <- DataFrame(sample_id = names(images),
                                           patient_id = patient_id,
                                           Response = Response)


#Cytomapper
cytomapper_IMC <- measureObjects(masks, image = images, img_id = "sample_id")

cytomapper_IMC



#Image visualization
library(SpatialExperiment)
library(cytomapper)

IMC <- readRDS('/Users/../IMC.rds')
images <- readRDS('/Users/../images.rds')
masks <- readRDS('/Users/../masks.rds')

# load filtered objects
IMC_filtered <- readRDS('/Users/../IMC_Filtered_UK_removed.rds')

# Ensure IMC_filtered$sample_id exists
# Assuming IMC_filtered is derived from IMC
IMC_filtered <- IMC  # Replace with your filtering logic if applicable
unique_sample_ids <- unique(IMC_filtered$sample_id)

# Filter the images and masks based on unique_sample_ids
images_filtered <- images[names(images) %in% unique_sample_ids]
masks_filtered <- masks[names(masks) %in% unique_sample_ids]

# Check the filtered results
images_filtered
masks_filtered


# Sample images
set.seed(220517)
# Predefined IDs you want to include

# Remaining IDs (excluding the specific ones)
res_id <- setdiff(unique(IMC$sample_id), non_res_id)

# For Non Responder Patients
#cur_id <- sample(unique(IMC$sample_id), 9)
cur_images <- images[names(images) %in% non_res_id] # non_responder patients
cur_masks <- masks[names(masks) %in% non_res_id] # non_responder patients

# For Responder Patients
cur_images <- images[names(images) %in% res_id] # responder patients
cur_masks <- masks[names(masks) %in% res_id] # responder patients

#Pixel visualization APOE+ TERM2
plotPixels(cur_images, 
           colour_by = c("S100","CD45", "CD31"),
           bcg = list(S100 = c(0, 5, 1),
                      CD45 = c(0, 5, 1),
                      CD31 = c(0, 5, 1)))

#Pixel visualization APOE+ TERM2
plotPixels(cur_images, 
           colour_by = c("CCL2","CCL4", "CD68"),
           bcg = list(CCL2 = c(0, 20, 1),
                      CCL4 = c(0, 2, 1),
                      CD68 = c(0, 10, 1)),
           colour = list(CCL2 = c("black", "white"),
                         CCL4 = c("black", "#FF0000"),
                         CD68 = c("black", "#00FF00")))


#Pixel visualization APOE+ TERM2
plotPixels(cur_images, 
           colour_by = c("CD8a", "TREM2","CD163", "APOE", "TIGIT"),
           bcg = list(TREM2 = c(0, 5, 1),
                      CD8a = c(0, 5, 1),
                      APOE = c(0, 1, 1),
                      TIGIT = c(0, 5, 1)),
           colour = list(CD8a = c("black", "#f79256"),
                         TREM2 = c("black", "#c1121f"),
                         CD163 = c("black", "#a7c957"),
                         APOE = c("black", "#00b4d8"),
                         TIGIT = c("black", "#f2f2f2")))


# ONLY TREM2
#Pixel visualization
plotPixels(cur_images, 
           colour_by = "HLADR",
           bcg = list(HLADR = c(0, 5, 1)),
           colour = list(HLADR = c("black", "cyan2")))


# visualize cells
plotCells(cur_masks,
          object = IMC, 
          cell_id = "ObjectNumber", 
          img_id = "sample_id",
          colour_by = c("CD8a", "TREM2", "CD163"),
          exprs_values = "exprs",
          colour = list(CD8a = c("black", "burlywood1"),
                        TREM2 = c("black", "cyan2"),
                        CD163 = c("black", "firebrick1")))


# define cell colors
# Get unique levels
celltype_levels <- unique(IMC_filtered$cluster_celltype)

# Define colors for each level
color_celltype_detail <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#ffff33", "#999999", "#003049", "#8dd3c7"
)

# Verify lengths match
if (length(color_celltype_detail) < length(celltype_levels)) {
  stop("Not enough colors for celltype levels.")
}

# Assign colors dynamically
names(color_celltype_detail) <- celltype_levels
colour_list <- list(cluster_celltype = color_celltype_detail)

# plot cells
plotCells(cur_masks,
          object = IMC_filtered, 
          cell_id = "ObjectNumber", 
          img_id = "sample_id",
          colour_by = "celltype_broad")


plotPixels(
  image = cur_images,
  mask = cur_masks,
  object = IMC_filtered,
  cell_id = "ObjectNumber",
  img_id = "sample_id",
  colour_by = c("S100", "CD45", "CD31"),
  outline_by = "celltype_braod",
  bcg = list(S100 = c(0, 5, 1), CD45 = c(0, 5, 1), CD31 = c(0, 5, 1)),
  colour = colour_list,  # Fixed color mapping
  thick = TRUE
)


############################################# UMAP ##########################################################################

# load IMC object if necessary
IMC_filtered <- readRDS('/Users/../IMC_Filtered_UK_removed.rds')

# load library
library(dittoSeq)
library(scater)
library(patchwork)
library(cowplot)
library(viridis)


# run UMAP
set.seed(220225)
IMC_filtered <- runUMAP(IMC_filtered, subset_row = rowData(IMC_filtered)$use_channel, exprs_values = "exprs") 

# run TSNE
IMC_filtered <- runTSNE(IMC_filtered, subset_row = rowData(IMC_filtered)$use_channel, exprs_values = "exprs") 

## UMAP colored by cell type and expression - dittoDimPlot
p1 <- dittoDimPlot(IMC_filtered, 
                   var = "patient_ID", 
                   reduction.use = "UMAP", 
                   size = 0.2,
                   do.label = FALSE) +
  scale_color_manual(values = metadata(IMC_filtered)$color_vectors$Patient_ID) +
  theme(legend.title = element_blank()) +
  ggtitle("Patient IDs")

p2 <- dittoDimPlot(IMC_filtered, 
                   var = "ROI", 
                   reduction.use = "UMAP", 
                   size = 0.2,
                   do.label = FALSE) +
  scale_color_manual(values = metadata(IMC_filtered)$color_vectors$ROI) +
  theme(legend.title = element_blank()) +
  ggtitle("ROIs")


p3 <- dittoDimPlot(IMC_filtered, 
                   var = "celltype_broad", 
                   reduction.use = "UMAP", 
                   size = 0.2,
                   do.label = FALSE) +
  scale_color_manual(values = metadata(IMC_filtered)$color_vectors$celltype_braod) +
  theme(legend.title = element_blank()) +
  ggtitle("Broad Celltype")

p4 <- dittoDimPlot(IMC_filtered, 
                   var = "cluster_celltype", 
                   reduction.use = "UMAP", 
                   size = 0.2,
                   do.label = FALSE) +
  scale_color_manual(values = metadata(IMC_filtered)$color_vectors$cluster_celltype) +
  theme(legend.title = element_blank()) +
  ggtitle("Detailed Celltype")

p1 + p2 + p3 + p4

# map individual genes
p5 <- dittoDimPlot(IMC_filtered, 
                   var = c("CD204", "CCL2", "CCL4"),  
                   assay = "exprs",
                   reduction.use = "UMAP", 
                   size = 0.2, 
                   colors = viridis(100), 
                   do.label = FALSE) +
  scale_color_viridis()

p5

################################### Subsetting Macrophage #############################################

# load IMC object if necessary
IMC_filtered <- readRDS('/Users/../IMC_Filtered.rds')

# subset based on response column 
Macrophage <- IMC_filtered[, IMC_filtered$celltype_broad == "Macrophage"] 

# define colors
color_vectors <- list()

my_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
               "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
               "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", "#393b79",
               "#637939", "#8c6d31", "#843c39", "#7b4173", "#3182bd", "#6baed6", "#9e9ac8",
               "#fd8d3c", "#e6550d", "#31a354", "#756bb1", "#636363", "#e41a1c", "#377eb8",
               "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999",
               "#fb8072", "#80b1d3", "#b3de69")

ROI <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
         "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
         "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", "#393b79",
         "#637939", "#8c6d31", "#843c39", "#7b4173", "#3182bd", "#6baed6", "#9e9ac8",
         "#fd8d3c", "#e6550d", "#31a354", "#756bb1", "#636363", "#e41a1c", "#377eb8",
         "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999",
         "#fb8072", "#80b1d3", "#b3de69")

sample <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
            "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
            "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", "#393b79",
            "#637939", "#8c6d31", "#843c39", "#7b4173", "#3182bd", "#6baed6", "#9e9ac8",
            "#fd8d3c", "#e6550d", "#31a354", "#756bb1", "#636363", "#e41a1c", "#377eb8",
            "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999",
            "#fb8072", "#80b1d3", "#b3de69")

Response <- c("#d6604d", "#4393c3")

Patient_ID <- c("#ff7f0e", "#2ca02c", "#d62728", "#ffff33", "#f781bf", "#3182bd")

clusters <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
              "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
              "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", "#393b79",
              "#637939", "#8c6d31", "#843c39", "#7b4173", "#3182bd", "#9e9ac8")

# add colors to metadata
color_vectors$ROI <- ROI
color_vectors$sample <- sample
color_vectors$Response <- Response
color_vectors$Patient_ID <- Patient_ID
color_vectors$clusters <- clusters
#color_vectors$immunetype <- immuntype
metadata(Macrophage)$color_vectors <- color_vectors


genes_for_heatmap <- c(
  "CD68", "HLADR", "CD15", "CD14", "CD163", "CD11b", 
  "CCL4", "CCL2", "TREM2", "APOE", "CD11c", "CD204", "CD38"
)

# Update rowData(IMC)$use_channel_heat
rowData(Macrophage)$use_channel <- rownames(Macrophage) %in% genes_for_heatmap

# select cells
cur_cells <- sample(seq_len(ncol(Macrophage)), 20000)

# get expression matrix
mat <- t(assay(Macrophage, "exprs")[rowData(Macrophage)$use_channel,])

# #scale rowwise
# scaled_mat <- t(scale(t(mat)))
# any(is.na(scaled_mat))   # Check for NA
# any(is.nan(scaled_mat))  # Check for NaN
# any(is.infinite(scaled_mat))  # Check for Inf
# which(is.na(scaled_mat), arr.ind = TRUE)   # Location of NAs
# which(is.nan(scaled_mat), arr.ind = TRUE)  # Location of NaNs
# scaled_mat[is.na(scaled_mat)] <- 0
# scaled_mat[is.nan(scaled_mat)] <- 0

# run phenograph
set.seed(230619)
out <- Rphenograph(mat, k = 45)

# Convert the clustering results into a vector
clusters <- factor(membership(out[[2]]))

# add a new column in IMC object 
Macrophage$clusters <- clusters

# Define a custom color palette
custom_colors <- colorRampPalette(c("blue", "white", "red"))(50)

# Set breaks for the heatmap to range from -2 to +2
breaks_seq <- seq(-2, 2, length.out = 50)


hm <- dittoHeatmap(Macrophage[,cur_cells], 
                   genes = rownames(Macrophage) [rowData(Macrophage)$use_channel],
                   assay = "exprs", scale = "row",
                   heatmap.colors = custom_colors,  # Custom color palette
                   breaks = breaks_seq,  # Breaks from -2 to +2
                   annot.by = c("clusters", "patient_ID", "Response"),
                   annot.colors = c(dittoColors(1)[1:length(unique(Macrophage$clusters))],
                                    metadata(Macrophage)$color_vectors$sample,
                                    metadata(Macrophage)$color_vectors$Response,
                                    metadata(Macrophage)$color_vectors$Patient_ID))


hm

# plot heatmap only marker genes
library(scuttle)

# Use aggregateAcrossCells with proper metadata propagation
celltype_mean <- aggregateAcrossCells(
  as(Macrophage, "SingleCellExperiment"), 
  ids = Macrophage$clusters,  # Group by cell type
  statistics = "mean",
  use.assay.type = "exprs",
  subset.row = genes_for_heatmap
)

# Create the heatmap
dittoHeatmap(
  celltype_mean,
  assay = "exprs", 
  cluster_cols = TRUE, 
  scale = "row",
  scaled.to.max = FALSE,
  heatmap.colors = custom_colors,  # Custom color palette
  annot.by = "cluster_celltype"  # Ensure these exist in colDat# Corrected annotation colors
)



############################################## UMAP ######################################################################

# run UMAP
set.seed(220225)
Macrophage <- runUMAP(Macrophage, subset_row = rowData(Macrophage)$use_channel, exprs_values = "exprs") 

# run TSNE
Macrophage <- runTSNE(Macrophage, subset_row = rowData(Macrophage)$use_channel, exprs_values = "exprs") 

## UMAP colored by cell type and expression - dittoDimPlot
p1 <- dittoDimPlot(Macrophage, 
                   var = "patient_ID", 
                   reduction.use = "UMAP", 
                   size = 0.2,
                   do.label = FALSE) +
  scale_color_manual(values = metadata(Macrophage)$color_vectors$Patient_ID) +
  theme(legend.title = element_blank()) +
  ggtitle("Patient IDs")

p2 <- dittoDimPlot(Macrophage, 
                   var = "ROI", 
                   reduction.use = "UMAP", 
                   size = 0.2,
                   do.label = FALSE) +
  scale_color_manual(values = metadata(Macrophage)$color_vectors$ROI) +
  theme(legend.title = element_blank()) +
  ggtitle("ROIs")


p3 <- dittoDimPlot(Macrophage, 
                   var = "clusters", 
                   reduction.use = "UMAP", 
                   size = 0.2,
                   do.label = TRUE) +
  scale_color_manual(values = metadata(Macrophage)$color_vectors$clusters) +
  theme(legend.title = element_blank()) +
  ggtitle("Myeloid Cluster")

p1 + p2 + p3 

# map individual genes
p5 <- dittoDimPlot(Macrophage, 
                   var = c("CD204", "CCL2", "CCL4"),  
                   assay = "exprs",
                   reduction.use = "UMAP", 
                   size = 0.2, 
                   colors = viridis(100), 
                   do.label = FALSE) +
  scale_color_viridis()

p5

# remove dendritic cells
# Define clusters to remove
clusters_to_remove <- c(5, 23, 25, 26, 27, 13, 24)

# Subset the Macrophage object to exclude these clusters
Macrophage_filtered <- Macrophage[, !Macrophage$clusters %in% clusters_to_remove]

# Check the unique clusters in the filtered object
unique(Macrophage_filtered$clusters)


###################################### cell type proportion in responder and non-responder patients #################################################
# load IMC object if necessary
Macrophage <- readRDS("/Users/../macrophage.rds")

# subset based on response column 
Macrophage_non_responder <- Macrophage_filtered[, Macrophage_filtered$Response == "Non-Responder"] 
Macrophage_responder <- Macrophage_filtered[, Macrophage_filtered$Response == "Responder"]

# by sample_id - percentage nonresponder
non_res <- dittoBarPlot(Macrophage_non_responder, 
                        var = "clusters", 
                        group.by = "sample_id") +
  scale_fill_manual(values = metadata(Macrophage_non_responder)$color_vectors$clusters)

# by sample_id - percentage responder
res <- dittoBarPlot(Macrophage_responder, 
                    var = "clusters", 
                    group.by = "sample_id") +
  scale_fill_manual(values = metadata(Macrophage_responder)$color_vectors$clusters)

########################################## agreegate bar plots ##########################################
library(dplyr)
library(ggplot2)

# Extract relevant data for aggregation
non_responder_data <- as.data.frame(colData(Macrophage_non_responder))  # Extract column data
responder_data <- as.data.frame(colData(Macrophage_responder))          # Extract column data

# Add a "Response" column to both datasets
non_responder_data <- non_responder_data %>%
  mutate(Response = "Non-Responder")
responder_data <- responder_data %>%
  mutate(Response = "Responder")

# Combine the two datasets
combined_data <- bind_rows(non_responder_data, responder_data)

# Calculate cell type proportions within each Response group
stacked_data <- combined_data %>%
  group_by(Response, clusters) %>%
  summarize(cell_count = n(), .groups = "drop") %>%
  group_by(Response) %>%
  mutate(average_percentage = cell_count / sum(cell_count) * 100)

# Create the stacked bar plot
stacked_bar_plot <- ggplot(stacked_data, aes(x = Response, y = average_percentage, fill = clusters)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = metadata(Macrophage)$color_vectors$clusters) +
  theme_minimal() +
  labs(title = "Cell Type Proportions",
       x = "Response Category", y = "Average Percentage", fill = "clusters") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display the plot
print(stacked_bar_plot)

######################################### t-test and zitter plot #########################################
library(dplyr)
library(ggplot2)

# Extract relevant data
data <- as.data.frame(colData(IMC_filtered))

# Calculate percentages for each sample and cell type
sample_data <- data %>%
  group_by(sample_id, Response, cluster_celltype) %>%
  summarize(cell_count = n(), .groups = "drop") %>%
  group_by(sample_id, Response) %>%
  mutate(percentage = cell_count / sum(cell_count) * 100)

# Filter data for Macrophage_ALA and Macrophage_IFN
filtered_data <- sample_data %>%
  filter(cluster_celltype %in% c("Macrophage_ALA", "Macrophage_IFN"))

# Perform t-tests for Macrophage_ALA and Macrophage_IFN
t_test_ala <- t.test(percentage ~ Response, data = filtered_data %>% filter(cluster_celltype == "Macrophage_ALA"))
t_test_ifn <- t.test(percentage ~ Response, data = filtered_data %>% filter(cluster_celltype == "Macrophage_IFN"))

# Print t-test results
cat("T-Test for Macrophage_ALA:\n")
print(t_test_ala)
cat("\nT-Test for Macrophage_IFN:\n")
print(t_test_ifn)

# Create jittered box plot
jitter_box_plot <- ggplot(filtered_data, aes(x = Response, y = percentage, fill = cluster_celltype)) +
  geom_boxplot(alpha = 0.5) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, aes(color = cluster_celltype)) +
  facet_wrap(~cluster_celltype) +
  scale_fill_manual(values = c("Macrophage_ALA" = "blue", "Macrophage_IFN" = "red")) +
  scale_color_manual(values = c("Macrophage_ALA" = "blue", "Macrophage_IFN" = "red")) +
  theme_minimal() +
  labs(title = "Macrophage_ALA and Macrophage_IFN Percentages by Response",
       x = "Response Category", y = "Percentage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  annotate("text", x = 1.5, y = max(filtered_data$percentage) + 5,
           label = sprintf("ALA p=%.3f\nIFN p=%.3f", t_test_ala$p.value, t_test_ifn$p.value),
           color = "black", size = 5, hjust = 0.5)

# Display the plot
print(jitter_box_plot)

############################################################################################################################
############################################################################################################################
#################################### IRA K15 CLUSTERING ###################################################################
############################################################################################################################
############################################################################################################################

# load RData from IRA
load('/Users/../final.RData')

#add meta data
# Add additional metadata
meta <- read_xlsx('/Users/../metadata.xlsx')
meta <- meta[, -1]


# # We can set the colnames of the object to generate unique identifiers per cell:
colnames(Melanoma_filtered) <- paste0(Melanoma_filtered$sample_id, "_", Melanoma_filtered$ObjectNumber)
# 
# # add meta data to spatial object 
# # # It is also often the case that sample-CCC_TMAcific metadata are available externally.
# IMC$sample_id = as.character(meta$Patient_ID)
# # It is also often the case that sample-CCC_TMAcific metadata are available externally.
# meta2 <- summary(as.factor(IMC$sample_id)) %>% as.data.frame()
# meta2 <- meta2 %>% dplyr::mutate(ID = rownames(meta2))
# meta2 <- meta2 %>% dplyr::select(ID)

# load R object from Ira
#IMC$indication <- meta$Indication[match(IMC$sample_id, meta$sample_id)]
#CCC_TMA$sample <- meta$Sample_ID[match(CCC_TMA$sample_id, meta$sample_id)]
Melanoma_filtered$ROI <- meta$Tissue_ID[match(Melanoma_filtered$sample_id, meta$sample_id)]
Melanoma_filtered$Response <- meta$Response[match(Melanoma_filtered$sample_id, meta$sample_id)]
Melanoma_filtered$patient_ID <- meta$PIZ[match(Melanoma_filtered$sample_id, meta$sample_id)]

# define and add colors
# Example of using this palette in a plot

color_vectors <- list()

my_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
               "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
               "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", "#393b79",
               "#637939", "#8c6d31", "#843c39", "#7b4173", "#3182bd", "#6baed6", "#9e9ac8",
               "#fd8d3c", "#e6550d", "#31a354", "#756bb1", "#636363", "#e41a1c", "#377eb8",
               "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999",
               "#fb8072", "#80b1d3", "#b3de69")

ROI <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
         "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
         "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", "#393b79",
         "#637939", "#8c6d31", "#843c39", "#7b4173", "#3182bd", "#6baed6", "#9e9ac8",
         "#fd8d3c", "#e6550d", "#31a354", "#756bb1", "#636363", "#e41a1c", "#377eb8",
         "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999",
         "#fb8072", "#80b1d3", "#b3de69")

sample <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
            "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
            "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", "#393b79",
            "#637939", "#8c6d31", "#843c39", "#7b4173", "#3182bd", "#6baed6", "#9e9ac8",
            "#fd8d3c", "#e6550d", "#31a354", "#756bb1", "#636363", "#e41a1c", "#377eb8",
            "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999",
            "#fb8072", "#80b1d3", "#b3de69")

Response <- c("#d6604d", "#4393c3")

Patient_ID <- c("#ff7f0e", "#2ca02c", "#d62728", "#ffff33", "#f781bf", "#3182bd")

cluster_celltype_filtered <- c("#ff7f0e", "#2ca02c", "#d62728", "#ffff33", "#f781bf", "#3182bd",
                               "#843c39", "#fb8072", "#b3de69", "#393b79", "#8c6d31", "#999999",  "#bcbd22",
                               "#984ea3", "#9e9ac8", "#bcbd22", "#17becf","#a65628")


celltype_sc <- c("#ff7f0e", "#2ca02c", "#d62728", "#ffff33", "#f781bf", "#3182bd",
                 "#843c39", "#fb8072", "#b3de69", "#393b79", "#8c6d31", "#999999",  
                 "#bcbd22", "#984ea3", "#17becf")


color_vectors$ROI <- ROI
color_vectors$sample <- sample
color_vectors$Response <- Response
color_vectors$Patient_ID <- Patient_ID
color_vectors$cluster_celltype_filtered <- cluster_celltype_filtered
color_vectors$celltype_sc <- celltype_sc
metadata(Melanoma_filtered)$color_vectors <- color_vectors

#for heat maps
rowData(Melanoma_filtered)$use_channel_heat <- !grepl("DNA1|DNA2|ki67|Hoechst|HH3|BCatenin|cCasp3|CytC", rownames(Melanoma_filtered))
#for umap
rowData(Melanoma_filtered)$use_channel <- !grepl("DNA1|DNA2|ki67|Hoechst|HH3|BCatenin|cCasp3|CytC", rownames(Melanoma_filtered))

cur_cells <- sample(seq_len(ncol(Melanoma_filtered)), 20000)

# Define a custom color palette
custom_colors <- colorRampPalette(c("blue", "white", "red"))(50)

# Set breaks for the heatmap to range from -2 to +2
breaks_seq <- seq(-2, 2, length.out = 50)

hm <- dittoHeatmap(Melanoma_filtered[,cur_cells], 
                   genes = rownames(Melanoma_filtered) [rowData(Melanoma_filtered)$use_channel],
                   assay = "exprs", scale = "row",
                   heatmap.colors = custom_colors,  # Custom color palette
                   breaks = breaks_seq,  # Breaks from -2 to +2
                   annot.by = c("cluster_celltype_filtered", "patient_ID", "Response"),
                   annot.colors = c(dittoColors(1)[1:length(unique(Melanoma_filtered$cluster_celltype_filtered))],
                                    metadata(Melanoma_filtered)$color_vectors$sample,
                                    metadata(Melanoma_filtered)$color_vectors$Response,
                                    metadata(Melanoma_filtered)$color_vectors$Patient_ID))


hm

# Select genes based on the updated marker_class
genes_detailed_celltype <- c("CD11c", "CD3", "CD8a", "CD45RO", "CD45", "HLADR", "CD68", "CD14", "CD163", "CD204", "Vimentin", "PanKer",
                             "ECadherin", "SMA", "CD31", "CD34", "CCL2", "CCL4", "Collagen1", "APOE", "CD15", "CD20", "S100", "GranzymeB")

####################################### celltype mean aggregate Clusters ########################################################

# plot heatmap only marker genes
library(scuttle)

# Use aggregateAcrossCells with proper metadata propagation
celltype_mean <- aggregateAcrossCells(
  as(Melanoma_filtered, "SingleCellExperiment"), 
  ids = Melanoma_filtered$cluster_celltype_filtered,  # Group by cell type
  statistics = "mean",
  use.assay.type = "exprs", 
  subset.row = genes_detailed_celltype
)

# Create the heatmap
dittoHeatmap(
  celltype_mean,
  assay = "exprs", 
  cluster_cols = TRUE, 
  scale = "row",
  heatmap.colors = custom_colors,  # Custom color palette
  annot.by = "cluster_celltype_filtered"  # Ensure these exist in colDat# Corrected annotation colors
)

################################ final naming of cell types ########################################

# now create a broader cell type category
celltype_sc <- recode(Melanoma_filtered$cluster_celltype_filtered,
                      "Dendritic cells" = "Dendritic cells",
                      "T cells" = "T cells",
                      "CD8 T cells" = "CD8 T cells",
                      "Macrophages" = "Macrophages",
                      "Stroma" = "Stroma",
                      "Endothelial progenitor cells" = "Endothelial progenitor",
                      "CCL2 CCL4" = "CAF",
                      "Tumor cells" = "Melanoma",
                      "APOE cells" = "Keratinocytes",
                      "Granulocytes" ="Granulocytes",
                      "HSPCs" = "Undefined",
                      "CCL2 CCL4 HSPCs" = "Perivascular CAF",
                      "Undefined" = "Undefined",
                      "HSCPs" = "Undefined",
                      "APC" = "APC",
                      "B cells" = "B cells",
                      "EMT cells" = "EMT cells",
                      "CCL2 CCL4 APOE" = "CAF"
)

Melanoma_filtered$celltype_sc <- celltype_sc


####################### heatmaps with final cell types ##############################
#hm
hm <- dittoHeatmap(Melanoma_filtered[,cur_cells], 
                   genes = rownames(Melanoma_filtered) [rowData(Melanoma_filtered)$use_channel],
                   assay = "exprs", scale = "row",
                   heatmap.colors = custom_colors,  # Custom color palette
                   breaks = breaks_seq,  # Breaks from -2 to +2
                   annot.by = c("celltype_sc", "patient_ID", "Response"),
                   annot.colors = c(dittoColors(1)[1:length(unique(Melanoma_filtered$celltype_sc))],
                                    metadata(Melanoma_filtered)$color_vectors$sample,
                                    metadata(Melanoma_filtered)$color_vectors$Response,
                                    metadata(Melanoma_filtered)$color_vectors$Patient_ID))


hm

# aggregate
# Use aggregateAcrossCells with proper metadata propagation
celltype_mean <- aggregateAcrossCells(
  as(Melanoma_filtered, "SingleCellExperiment"), 
  ids = Melanoma_filtered$celltype_sc,  # Group by cell type
  statistics = "mean",
  use.assay.type = "exprs", 
  subset.row = genes_detailed_celltype
)

# Create the heatmap
dittoHeatmap(
  celltype_mean,
  assay = "exprs", 
  cluster_cols = TRUE, 
  scale = "row",
  heatmap.colors = custom_colors,  # Custom color palette
  annot.by = "celltype_sc"  # Ensure these exist in colDat# Corrected annotation colors
)


# remove undefined cells 
# Subset the Melanoma_filtered object to exclude "Undefined" cells
Melanoma_wo_ud <- Melanoma_filtered[, Melanoma_filtered$celltype_sc != "Undefined"]

# Verify that "Undefined" cells have been removed
unique(Melanoma_wo_ud$celltype_sc)



################################## UMAP ###########################################################
# load IMC object if necessary
Melanoma_wo_ud <- readRDS('/Users/../melanoma_filtered_wo_ud.rds')

# load library
library(dittoSeq)
library(scater)
library(patchwork)
library(cowplot)
library(viridis)

#for heat maps
#rowData(Melanoma_filtered)$use_channel_heat <- !grepl("DNA1|DNA2|ki67|Hoechst|HH3|BCatenin|cCasp3|CytC", rownames(Melanoma_filtered))
#for umap
rowData(Melanoma_wo_ud)$use_channel <- !grepl("DNA1|DNA2|ki67|Hoechst|HH3|BCatenin|cCasp3|CytC", rownames(Melanoma_wo_ud))

# run UMAP
set.seed(220225)
Melanoma_wo_ud <- runUMAP(Melanoma_wo_ud, subset_row = rowData(Melanoma_wo_ud)$use_channel, exprs_values = "exprs", scale= TRUE) 

# # run TSNE
# IMC_filtered <- runTSNE(IMC_filtered, subset_row = rowData(IMC_filtered)$use_channel, exprs_values = "exprs") 

## UMAP colored by cell type and expression - dittoDimPlot
p1 <- dittoDimPlot(Melanoma_wo_ud, 
                   var = "patient_ID", 
                   reduction.use = "UMAP", 
                   size = 0.2,
                   do.label = FALSE) +
  scale_color_manual(values = metadata(Melanoma_wo_ud)$color_vectors$Patient_ID) +
  theme(legend.title = element_blank()) +
  ggtitle("Patient IDs")

p2 <- dittoDimPlot(Melanoma_wo_ud, 
                   var = "ROI", 
                   reduction.use = "UMAP", 
                   size = 0.2,
                   do.label = FALSE) +
  scale_color_manual(values = metadata(Melanoma_wo_ud)$color_vectors$ROI) +
  theme(legend.title = element_blank()) +
  ggtitle("ROIs")


p3 <- dittoDimPlot(Melanoma_wo_ud, 
                   var = "celltype_sc", 
                   reduction.use = "UMAP", 
                   size = 0.2,
                   do.label = FALSE) +
  scale_color_manual(values = metadata(Melanoma_wo_ud)$color_vectors$celltype_sc) +
  theme(legend.title = element_blank()) +
  ggtitle("CelltypeSC")

p4 <- dittoDimPlot(Melanoma_wo_ud, 
                   var = "Response", 
                   reduction.use = "UMAP", 
                   size = 0.2,
                   do.label = FALSE) +
  scale_color_manual(values = metadata(Melanoma_wo_ud)$color_vectors$Response) +
  theme(legend.title = element_blank()) +
  ggtitle("Response")

p1 + p2 + p3 + p4

# map individual genes
p5 <- dittoDimPlot(Melanoma_wo_ud, 
                   var = c("CD204", "CCL2", "CCL4"),  
                   assay = "exprs",
                   reduction.use = "UMAP", 
                   size = 0.2, 
                   colors = viridis(100), 
                   do.label = FALSE) +
  scale_color_viridis()

p5

# Define the genes to plot
genes_to_plot <- c("CD3", "CD8a", "CD20", "CD204", "CD163", "CD68", "CD11c", "SMA", "ECadherin", "S100", "CD31")

# Check which genes are present in the object
missing_genes <- setdiff(genes_to_plot, rownames(Melanoma_wo_ud))
present_genes <- intersect(genes_to_plot, rownames(Melanoma_wo_ud))

# Use only the genes that are present
genes_to_plot <- present_genes

# Generate UMAP plots
plot_list <- lapply(genes_to_plot, function(gene) {
  plotReducedDim(
    Melanoma_wo_ud,
    dimred = "UMAP",
    colour_by = gene,
    by_exprs_values = "exprs",
    point_size = 0.2
  )
})

# Combine all plots into a grid
plot_grid(plotlist = plot_list)

# violin plot
# Violin Plot - plotExpression
vp <- plotExpression(
  Melanoma_wo_ud,
  features = genes_to_plot,  # Use the specified genes
  x = "celltype_sc",           # Group by cell type
  exprs_values = "exprs",   # Expression values to plot
  colour_by = "celltype_sc"    # Color by cell type
) + 
  theme(axis.text.x = element_text(angle = 90)) +  # Rotate x-axis labels
  scale_color_manual(values = metadata(Melanoma_wo_ud)$color_vectors$celltype)

##################################### MACROPHAGE SUBCLUSTER SPLITTING ###############################################

################################### Subsetting Macrophage #############################################

# load IMC object if necessary
Melanoma_wo_ud <- readRDS('/Users/../melanoma_filtered_wo_ud.rds')

# subset based on response column 
Macrophage <- Melanoma_wo_ud[, Melanoma_wo_ud$celltype_sc == "Macrophages"] 

# define colors
color_vectors <- list()

my_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
               "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
               "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", "#393b79",
               "#637939", "#8c6d31", "#843c39", "#7b4173", "#3182bd", "#6baed6", "#9e9ac8",
               "#fd8d3c", "#e6550d", "#31a354", "#756bb1", "#636363", "#e41a1c", "#377eb8",
               "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999",
               "#fb8072", "#80b1d3", "#b3de69")

ROI <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
         "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
         "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", "#393b79",
         "#637939", "#8c6d31", "#843c39", "#7b4173", "#3182bd", "#6baed6", "#9e9ac8",
         "#fd8d3c", "#e6550d", "#31a354", "#756bb1", "#636363", "#e41a1c", "#377eb8",
         "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999",
         "#fb8072", "#80b1d3", "#b3de69")

sample <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
            "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
            "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", "#393b79",
            "#637939", "#8c6d31", "#843c39", "#7b4173", "#3182bd", "#6baed6", "#9e9ac8",
            "#fd8d3c", "#e6550d", "#31a354", "#756bb1", "#636363", "#e41a1c", "#377eb8",
            "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999",
            "#fb8072", "#80b1d3", "#b3de69")

Response <- c("#d6604d", "#4393c3")

Patient_ID <- c("#ff7f0e", "#2ca02c", "#d62728", "#ffff33", "#f781bf", "#3182bd")

clusters <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
              "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
              "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", "#393b79",
              "#637939", "#8c6d31", "#843c39", "#7b4173", "#3182bd", "#9e9ac8")


cluster_celltype <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
                      "#7f7f7f", "#bcbd22", "#17becf", "#f7b6d2")


# add colors to metadata
color_vectors$ROI <- ROI
color_vectors$sample <- sample
color_vectors$Response <- Response
color_vectors$Patient_ID <- Patient_ID
color_vectors$clusters <- clusters
color_vectors$cluster_celltype <- cluster_celltype
#color_vectors$immunetype <- immuntype
metadata(Macrophage_filtered)$color_vectors <- color_vectors

genes_for_heatmap <- c(
  "CD68", "HLADR", "CD15", "CD14", "CD163", "CD11b", 
  "CCL4", "CCL2", "TREM2", "APOE", "CD11c", "CD204", "CD38"
)

# Update rowData(IMC)$use_channel_heat
rowData(Macrophage)$use_channel <- rownames(Macrophage) %in% genes_for_heatmap

# select cells
cur_cells <- sample(seq_len(ncol(Macrophage)), 20000)

# get expression matrix
mat <- t(assay(Macrophage, "exprs")[rowData(Macrophage)$use_channel,])

# #scale rowwise
# scaled_mat <- t(scale(t(mat)))
# any(is.na(scaled_mat))   # Check for NA
# any(is.nan(scaled_mat))  # Check for NaN
# any(is.infinite(scaled_mat))  # Check for Inf
# which(is.na(scaled_mat), arr.ind = TRUE)   # Location of NAs
# which(is.nan(scaled_mat), arr.ind = TRUE)  # Location of NaNs
# scaled_mat[is.na(scaled_mat)] <- 0
# scaled_mat[is.nan(scaled_mat)] <- 0

# run phenograph
set.seed(230619)
out <- Rphenograph(mat, k = 205)

# Convert the clustering results into a vector
clusters <- factor(membership(out[[2]]))

# add a new column in IMC object 
Macrophage$clusters <- clusters

# Define a custom color palette
custom_colors <- colorRampPalette(c("blue", "white", "red"))(50)

# Set breaks for the heatmap to range from -2 to +2
breaks_seq <- seq(-2, 2, length.out = 50)


hm <- dittoHeatmap(Macrophage[,cur_cells], 
                   genes = rownames(Macrophage) [rowData(Macrophage)$use_channel],
                   assay = "exprs", scale = "row",
                   heatmap.colors = custom_colors,  # Custom color palette
                   breaks = breaks_seq,  # Breaks from -2 to +2
                   annot.by = c("clusters", "patient_ID", "Response"),
                   annot.colors = c(dittoColors(1)[1:length(unique(Macrophage$clusters))],
                                    metadata(Macrophage)$color_vectors$sample,
                                    metadata(Macrophage)$color_vectors$Response,
                                    metadata(Macrophage)$color_vectors$Patient_ID))


hm

# plot heatmap only marker genes
library(scuttle)

# Use aggregateAcrossCells with proper metadata propagation
celltype_mean <- aggregateAcrossCells(
  as(Macrophage, "SingleCellExperiment"), 
  ids = Macrophage$clusters,  # Group by cell type
  statistics = "mean",
  use.assay.type = "exprs",
  subset.row = genes_for_heatmap
)

# Create the heatmap
dittoHeatmap(
  celltype_mean,
  assay = "exprs", 
  cluster_cols = TRUE, 
  scale = "row",
  scaled.to.max = FALSE,
  heatmap.colors = custom_colors,  # Custom color palette
  annot.by = "clusters"  # Ensure these exist in colDat# Corrected annotation colors
)

# annotate clusters
library(dplyr)
cluster_celltype <- recode(Macrophage$clusters,
                           "1" = "Monocytes",
                           "2" = "Myeloid",
                           "3" = "LATAM",
                           "4" = "TAM",
                           "5" = "TAM",
                           "6" = "Antigen-Presenting Macrophages",
                           "7" = "APOE+ cell",
                           "8" = "Antigen-Presenting Macrophages",
                           "9" = "APC",
                           "10" = "Undefined",
                           "11" = "IFNTAM",
                           "12" = "Macrophages",
                           "13" = "Dendritic cells",
                           "14" = "Undefined",
                           '15' = "Undefined"
)

Macrophage$cluster_celltype <- cluster_celltype



############################################## UMAP ######################################################################

# run UMAP
set.seed(220225)
Macrophage_filtered <- runUMAP(Macrophage_filtered, subset_row = rowData(Macrophage_filtered)$use_channel, exprs_values = "exprs", scale = TRUE, n_neighbors = 10) 

# # run TSNE
# Macrophage <- runTSNE(Macrophage, subset_row = rowData(Macrophage)$use_channel, exprs_values = "exprs") 

## UMAP colored by cell type and expression - dittoDimPlot
p1 <- dittoDimPlot(Macrophage_filtered, 
                   var = "patient_ID", 
                   reduction.use = "UMAP", 
                   size = 0.2,
                   do.label = FALSE) +
  scale_color_manual(values = metadata(Macrophage_filtered)$color_vectors$Patient_ID) +
  theme(legend.title = element_blank()) +
  ggtitle("Patient IDs")

p2 <- dittoDimPlot(Macrophage_filtered, 
                   var = "ROI", 
                   reduction.use = "UMAP", 
                   size = 0.2,
                   do.label = FALSE) +
  scale_color_manual(values = metadata(Macrophage_filtered)$color_vectors$ROI) +
  theme(legend.title = element_blank()) +
  ggtitle("ROIs")


p3 <- dittoDimPlot(Macrophage_filtered, 
                   var = "cluster_celltype", 
                   reduction.use = "UMAP", 
                   size = 0.2,
                   do.label = FALSE) +
  scale_color_manual(values = metadata(Macrophage_filtered)$color_vectors$cluster_celltype) +
  theme(legend.title = element_blank()) +
  ggtitle("Myeloid Cluster")

p1 + p2 + p3 

# map individual genes
p5 <- dittoDimPlot(Macrophage_filtered, 
                   var = c("APOE", "TREM2", "CD68", "CCL4"),  
                   assay = "exprs",
                   reduction.use = "UMAP", 
                   size = 0.2, 
                   colors = viridis(100), 
                   do.label = FALSE) +
  scale_color_viridis()

p5

###################################### cell type proportion in responder and non-responder patients #################################################
# subset based on response column 
Macrophage_non_responder <- Macrophage_filtered[, Macrophage_filtered$Response == "Non-Responder"] 
Macrophage_responder <- Macrophage_filtered[, Macrophage_filtered$Response == "Responder"]

# by sample_id - percentage nonresponder
non_res <- dittoBarPlot(Macrophage_non_responder, 
                        var = "cluster_celltype", 
                        group.by = "sample_id") +
  scale_fill_manual(values = metadata(Macrophage_non_responder)$color_vectors$cluster_celltype)

# by sample_id - percentage responder
res <- dittoBarPlot(Macrophage_responder, 
                    var = "cluster_celltype", 
                    group.by = "sample_id") +
  scale_fill_manual(values = metadata(Macrophage_responder)$color_vectors$cluster_celltype)

########################################## agreegate bar plots ##########################################
library(dplyr)
library(ggplot2)

# Extract relevant data for aggregation
non_responder_data <- as.data.frame(colData(Macrophage_non_responder))  # Extract column data
responder_data <- as.data.frame(colData(Macrophage_responder))          # Extract column data

# Add a "Response" column to both datasets
non_responder_data <- non_responder_data %>%
  mutate(Response = "Non-Responder")
responder_data <- responder_data %>%
  mutate(Response = "Responder")

# Combine the two datasets
combined_data <- bind_rows(non_responder_data, responder_data)

# Calculate cell type proportions within each Response group
stacked_data <- combined_data %>%
  group_by(Response, cluster_celltype) %>%
  summarize(cell_count = n(), .groups = "drop") %>%
  group_by(Response) %>%
  mutate(average_percentage = cell_count / sum(cell_count) * 100)

# Create the stacked bar plot
stacked_bar_plot <- ggplot(stacked_data, aes(x = Response, y = average_percentage, fill = cluster_celltype)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = metadata(Macrophage)$color_vectors$cluster_celltype) +
  theme_minimal() +
  labs(title = "Cell Type Proportions",
       x = "Response Category", y = "Average Percentage", fill = "clusters") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display the plot
print(stacked_bar_plot)

######################################### t-test and zitter plot #########################################
library(dplyr)
library(ggplot2)

# Extract relevant data
data <- as.data.frame(colData(Macrophage_filtered))

# Calculate percentages for each sample and cell type
sample_data <- data %>%
  group_by(sample_id, Response, cluster_celltype) %>%
  summarize(cell_count = n(), .groups = "drop") %>%
  group_by(sample_id, Response) %>%
  mutate(percentage = cell_count / sum(cell_count) * 100)

# Filter data for Macrophage_ALA and Macrophage_IFN
filtered_data <- sample_data %>%
  filter(cluster_celltype %in% c("LATAM", "IFNTAM"))

# Perform t-tests for Macrophage_ALA and Macrophage_IFN
t_test_ala <- t.test(percentage ~ Response, data = filtered_data %>% filter(cluster_celltype == "LATAM"))
t_test_ifn <- t.test(percentage ~ Response, data = filtered_data %>% filter(cluster_celltype == "IFNTAM"))

# Print t-test results
cat("T-Test for Macrophage_ALA:\n")
print(t_test_ala)
cat("\nT-Test for Macrophage_IFN:\n")
print(t_test_ifn)

# Create jittered box plot
jitter_box_plot <- ggplot(filtered_data, aes(x = Response, y = percentage, fill = cluster_celltype)) +
  geom_boxplot(alpha = 0.5) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7, aes(color = cluster_celltype)) +
  facet_wrap(~cluster_celltype) +
  scale_fill_manual(values = c("Macrophage_ALA" = "blue", "Macrophage_IFN" = "red")) +
  scale_color_manual(values = c("Macrophage_ALA" = "blue", "Macrophage_IFN" = "red")) +
  theme_minimal() +
  labs(title = "Macrophage_ALA and Macrophage_IFN Percentages by Response",
       x = "Response Category", y = "Percentage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  annotate("text", x = 1.5, y = max(filtered_data$percentage) + 5,
           label = sprintf("ALA p=%.3f\nIFN p=%.3f", t_test_ala$p.value, t_test_ifn$p.value),
           color = "black", size = 5, hjust = 0.5)

# Display the plot
print(jitter_box_plot)

################################### transfer macrophage annotation to melanoma_filtered object ######################################

# Step 0: Create a mapping from `Macrophage_filtered`
annotation_mapping <- data.frame(
  cell_id = colnames(Macrophage),                # Cell IDs from Macrophage_filtered
  cluster_celltype = Macrophage$cluster_celltype # Annotations from Macrophage_filtered
)

# Step 1: Convert `celltype_sc2` to a character vector
Melanoma_wo_ud$celltype_sc2 <- as.character(Melanoma_wo_ud$celltype_sc)

# Step 2: Update `celltype_sc2` with the new annotations
Melanoma_wo_ud$celltype_sc2[!is.na(matched_indices)] <- as.character(annotation_mapping$cluster_celltype[matched_indices[!is.na(matched_indices)]])

# Step 3: (Optional) Convert back to a factor if needed
Melanoma_wo_ud$celltype_sc2 <- factor(Melanoma_wo_ud$celltype_sc2)

# Verify the updated annotations
table(Melanoma_wo_ud$celltype_sc2)

# Subset the spatial object to exclude "Unknown"
melanoma_wo_ud_macro <- Melanoma_wo_ud[, Melanoma_wo_ud$celltype_sc2 != "Undefined"]

# define colors
color_vectors <- list()

celltype_sc2 <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#8c6d31", "#9467bd", 
                  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
                  "#393b79", "#637939", "#d62728", "#843c39", "#7b4173", 
                  "#3182bd", "#6baed6", "#9e9ac8", "#fd8d3c", "#e6550d", 
                  "#31a354")


# add colors to metadata
color_vectors$celltype_sc2 <- celltype_sc2
#color_vectors$immunetype <- immuntype
metadata(melanoma_wo_ud_macro)$color_vectors <- color_vectors


# check ump
p3 <- dittoDimPlot(melanoma_wo_ud_macro, 
                   var = "celltype_sc2", 
                   reduction.use = "UMAP", 
                   size = 0.2,
                   do.label = FALSE) +
  scale_color_manual(values = metadata(melanoma_wo_ud_macro)$color_vectors$celltype_sc2) +
  theme(legend.title = element_blank()) +
  ggtitle("cell Cluster")

p3


##################################### CD8 T cell SUBCLUSTER SPLITTING ###############################################

################################### Subsetting CD8 #############################################

# subset based on response column 
CD8 <- melanoma_wo_ud_macro[, melanoma_wo_ud_macro$celltype_sc == "CD8 T cells"] 

# define colors
color_vectors <- list()

my_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
               "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
               "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", "#393b79",
               "#637939", "#8c6d31", "#843c39", "#7b4173", "#3182bd", "#6baed6", "#9e9ac8",
               "#fd8d3c", "#e6550d", "#31a354", "#756bb1", "#636363", "#e41a1c", "#377eb8",
               "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999",
               "#fb8072", "#80b1d3", "#b3de69")

ROI <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
         "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
         "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", "#393b79",
         "#637939", "#8c6d31", "#843c39", "#7b4173", "#3182bd", "#6baed6", "#9e9ac8",
         "#fd8d3c", "#e6550d", "#31a354", "#756bb1", "#636363", "#e41a1c", "#377eb8",
         "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999",
         "#fb8072", "#80b1d3", "#b3de69")

sample <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
            "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
            "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", "#393b79",
            "#637939", "#8c6d31", "#843c39", "#7b4173", "#3182bd", "#6baed6", "#9e9ac8",
            "#fd8d3c", "#e6550d", "#31a354", "#756bb1", "#636363", "#e41a1c", "#377eb8",
            "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999",
            "#fb8072", "#80b1d3", "#b3de69")

Response <- c("#d6604d", "#4393c3")

Patient_ID <- c("#ff7f0e", "#2ca02c", "#d62728", "#ffff33", "#f781bf", "#3182bd")

clusters <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
              "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
              "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5", "#393b79",
              "#637939", "#8c6d31", "#843c39", "#7b4173", "#3182bd", "#9e9ac8")


cluster_celltype <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
                      "#7f7f7f", "#bcbd22", "#17becf", "#f7b6d2")


# add colors to metadata
color_vectors$ROI <- ROI
color_vectors$sample <- sample
color_vectors$Response <- Response
color_vectors$Patient_ID <- Patient_ID
color_vectors$clusters <- clusters
color_vectors$cluster_celltype <- cluster_celltype
#color_vectors$immunetype <- immuntype
metadata(CD8_subset)$color_vectors <- color_vectors

genes_for_heatmap <- c(
  "CD3", "Tim3",
  "CD8a", "PD1", "GranzymeB"
)

# Update rowData(IMC)$use_channel_heat
rowData(CD8_subset)$use_channel <- rownames(CD8_subset) %in% genes_for_heatmap

# select cells
cur_cells <- sample(seq_len(ncol(CD8_subset)))

# get expression matrix
mat <- t(assay(CD8_subset, "exprs")[rowData(CD8)$use_channel,])

# #scale rowwise
# scaled_mat <- t(scale(t(mat)))
# any(is.na(scaled_mat))   # Check for NA
# any(is.nan(scaled_mat))  # Check for NaN
# any(is.infinite(scaled_mat))  # Check for Inf
# which(is.na(scaled_mat), arr.ind = TRUE)   # Location of NAs
# which(is.nan(scaled_mat), arr.ind = TRUE)  # Location of NaNs
# scaled_mat[is.na(scaled_mat)] <- 0
# scaled_mat[is.nan(scaled_mat)] <- 0

# run phenograph
set.seed(230619)
out <- Rphenograph(mat, k = 250)

# Convert the clustering results into a vector
clusters <- factor(membership(out[[2]]))

# add a new column in IMC object 
CD8_subset$clusters <- clusters

# Define a custom color palette
custom_colors <- colorRampPalette(c("blue", "white", "red"))(50)

# Set breaks for the heatmap to range from -2 to +2
breaks_seq <- seq(-2, 2, length.out = 50)


hm <- dittoHeatmap(CD8_subset[,cur_cells], 
                   genes = rownames(CD8_subset) [rowData(CD8_subset)$use_channel],
                   assay = "exprs", scale = "row",
                   heatmap.colors = custom_colors,  # Custom color palette
                   breaks = breaks_seq,  # Breaks from -2 to +2
                   annot.by = c("clusters", "patient_ID", "Response"),
                   annot.colors = c(dittoColors(1)[1:length(unique(CD8_subset$clusters))],
                                    metadata(CD8_subset)$color_vectors$sample,
                                    metadata(CD8_subset)$color_vectors$Response,
                                    metadata(CD8_subset)$color_vectors$Patient_ID))


hm


# Subset CD8 object to keep only clusters 2, 3, and 5
CD8_subset <- CD8_subset[, CD8_subset$clusters %in% c(2, 3, 5)]

# Verify the retained clusters
unique(CD8_subset$clusters)


# plot heatmap only marker genes
library(scuttle)

# Use aggregateAcrossCells with proper metadata propagation
celltype_mean <- aggregateAcrossCells(
  as(CD8, "SingleCellExperiment"), 
  ids = CD8$clusters,  # Group by cell type
  statistics = "mean",
  use.assay.type = "exprs",
  subset.row = genes_for_heatmap
)

# Create the heatmap
dittoHeatmap(
  celltype_mean,
  assay = "exprs", 
  cluster_cols = TRUE, 
  scale = "row",
  scaled.to.max = FALSE,
  heatmap.colors = custom_colors,  # Custom color palette
  annot.by = "clusters"  # Ensure these exist in colDat# Corrected annotation colors
)

#load
CD8 <- readRDS("/Users/...CD8.rds")

# annotate clusters
library(dplyr)
cluster_celltype <- recode(CD8$clusters,
                           "1" = "TIM+ CD8",
                           "2" = "CD8 T cells",
                           "3" = "Gymb+ CD8",
                           "4" = "CD8 T cells",
                           "5" = "CD8 T cells",
                           "6" = "CD8 T cells",
                           "7" = "CD8 T cells",
                           "8" = "CD8 T cells"
)

CD8$cluster_celltype <- cluster_celltype

hm <- dittoHeatmap(CD8[,cur_cells], 
                   genes = rownames(CD8) [rowData(CD8)$use_channel],
                   assay = "exprs", scale = "row",
                   heatmap.colors = custom_colors,  # Custom color palette
                   breaks = breaks_seq,  # Breaks from -2 to +2
                   annot.by = c("cluster_celltype", "patient_ID", "Response"),
                   annot.colors = c(dittoColors(1)[1:length(unique(CD8$cluster_celltype))],
                                    metadata(CD8)$color_vectors$sample,
                                    metadata(CD8)$color_vectors$Response,
                                    metadata(CD8)$color_vectors$Patient_ID))


hm

# Violin Plot - plotExpression
plotExpression(CD8[,cur_cells], 
               features = genes_for_heatmap,
               x = "cluster_celltype", 
               exprs_values = "exprs", 
               colour_by = "cluster_celltype") +
  theme(axis.text.x =  element_text(angle = 90))+
  scale_color_manual(values = metadata(CD8)$color_vectors$cluster_celltype)


# Verify that "Undefined" cells have been removed
unique(CD8_filtered$cluster_celltype)

############################################## UMAP ######################################################################
# run UMAP
set.seed(220225)
CD8 <- runUMAP(CD8, subset_row = rowData(CD8)$use_channel, exprs_values = "exprs") 

# # run TSNE
# Macrophage <- runTSNE(Macrophage, subset_row = rowData(Macrophage)$use_channel, exprs_values = "exprs") 

## UMAP colored by cell type and expression - dittoDimPlot
p1 <- dittoDimPlot(CD8, 
                   var = "patient_ID", 
                   reduction.use = "UMAP", 
                   size = 0.2,
                   do.label = FALSE) +
  scale_color_manual(values = metadata(CD8)$color_vectors$Patient_ID) +
  theme(legend.title = element_blank()) +
  ggtitle("Patient IDs")

p2 <- dittoDimPlot(CD8, 
                   var = "ROI", 
                   reduction.use = "UMAP", 
                   size = 0.2,
                   do.label = FALSE) +
  scale_color_manual(values = metadata(CD8)$color_vectors$ROI) +
  theme(legend.title = element_blank()) +
  ggtitle("ROIs")


p3 <- dittoDimPlot(CD8, 
                   var = "cluster_celltype", 
                   reduction.use = "UMAP", 
                   size = 0.2,
                   do.label = FALSE) +
  scale_color_manual(values = metadata(CD8)$color_vectors$cluster_celltype) +
  theme(legend.title = element_blank()) +
  ggtitle("CD8 Cluster")

p4 <- dittoDimPlot(CD8, 
                   var = "Response", 
                   reduction.use = "UMAP", 
                   size = 0.2,
                   do.label = FALSE) +
  scale_color_manual(values = metadata(CD8)$color_vectors$Response) +
  theme(legend.title = element_blank()) +
  ggtitle("CD8 Cluster")



p1 + p2 + p3 + p4

# map individual genes
p5 <- dittoDimPlot(CD8_filtered, 
                   var = c("CD8a", "GranzymeB", "Tim3"),  
                   assay = "exprs",
                   reduction.use = "UMAP", 
                   size = 0.2, 
                   colors = viridis(100), 
                   do.label = FALSE) +
  scale_color_viridis()

p5

p6 <- dittoDimPlot(CD8_filtered, 
                   var = "Tim3",  
                   assay = "exprs",
                   reduction.use = "UMAP", 
                   size = 0.2, 
                   colors = viridis(100), 
                   do.label = FALSE) +
  scale_color_viridis()



###################################### cell type proportion in responder and non-responder patients #################################################
# subset based on response column 
CD8_non_responder <- CD8[, CD8$Response == "Non-Responder"] 
CD8_responder <- CD8[, CD8$Response == "Responder"]

# by sample_id - percentage nonresponder
non_res <- dittoBarPlot(CD8_non_responder, 
                        var = "cluster_celltype", 
                        group.by = "sample_id") +
  scale_fill_manual(values = metadata(CD8_non_responder)$color_vectors$cluster_celltype)

# by sample_id - percentage responder
res <- dittoBarPlot(CD8_responder, 
                    var = "cluster_celltype", 
                    group.by = "sample_id") +
  scale_fill_manual(values = metadata(CD8_responder)$color_vectors$cluster_celltype)

########################################## agreegate bar plots ##########################################
library(dplyr)
library(ggplot2)

# Extract relevant data for aggregation
non_responder_data <- as.data.frame(colData(CD8_non_responder))  # Extract column data
responder_data <- as.data.frame(colData(CD8_responder))          # Extract column data

# Add a "Response" column to both datasets
non_responder_data <- non_responder_data %>%
  mutate(Response = "Non-Responder")
responder_data <- responder_data %>%
  mutate(Response = "Responder")

# Combine the two datasets
combined_data <- bind_rows(non_responder_data, responder_data)

# Calculate cell type proportions within each Response group
stacked_data <- combined_data %>%
  group_by(Response, cluster_celltype) %>%
  summarize(cell_count = n(), .groups = "drop") %>%
  group_by(Response) %>%
  mutate(average_percentage = cell_count / sum(cell_count) * 100)

# Create the stacked bar plot
stacked_bar_plot <- ggplot(stacked_data, aes(x = Response, y = average_percentage, fill = cluster_celltype)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = metadata(CD8_filtered)$color_vectors$cluster_celltype) +
  theme_minimal() +
  labs(title = "Cell Type Proportions",
       x = "Response Category", y = "Average Percentage", fill = "clusters") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display the plot
print(stacked_bar_plot)

################################### transfer macrophage annotation to melanoma_filtered object ######################################

# load the base object to which the CD8 annotations will be transfered 

# check cell types # macrophage annotation should be present
unique(melanoma_wo_ud_macro$celltype_sc2)

# check the cell types of CD8 object 
unique(CD8$cluster_celltype)

# Step 0: Create a mapping from `CD8`
annotation_mapping <- data.frame(
  cell_id = colnames(CD8),                # Cell IDs from CD8
  cluster_celltype = CD8$cluster_celltype # Annotations from CD8
)

# Step 1: Create a new column `celltype_sc3` and convert it to character
melanoma_wo_ud_macro$celltype_sc3 <- as.character(melanoma_wo_ud_macro$celltype_sc2)

# Step 2: Find matched indices
matched_indices <- match(colnames(melanoma_wo_ud_macro), annotation_mapping$cell_id)

# Step 3: Update `celltype_sc3` with the new annotations
melanoma_wo_ud_macro$celltype_sc3[!is.na(matched_indices)] <- as.character(
  annotation_mapping$cluster_celltype[matched_indices[!is.na(matched_indices)]]
)

# Step 4: (Optional) Convert `celltype_sc3` back to a factor
melanoma_wo_ud_macro$celltype_sc3 <- factor(melanoma_wo_ud_macro$celltype_sc3)

# Verify the updated annotations
table(melanoma_wo_ud_macro$celltype_sc3)

# define colors
color_vectors <- list()

celltype_sc2 <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#8c6d31", "#9467bd", 
                  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
                  "#393b79", "#637939", "#d62728", "#843c39", "#7b4173", 
                  "#3182bd", "#6baed6", "#9e9ac8", "#fd8d3c", "#e6550d", 
                  "#31a354")


# add colors to metadata
color_vectors$celltype_sc2 <- celltype_sc2
#color_vectors$immunetype <- immuntype
metadata(melanoma_wo_ud_macro)$color_vectors <- color_vectors


# check ump
p3 <- dittoDimPlot(melanoma_wo_ud_macro, 
                   var = "celltype_sc2", 
                   reduction.use = "UMAP", 
                   size = 0.2,
                   do.label = FALSE) +
  scale_color_manual(values = metadata(melanoma_wo_ud_macro)$color_vectors$celltype_sc2) +
  theme(legend.title = element_blank()) +
  ggtitle("cell Cluster")

p3

##################################### MACROPHAGE SUBCLUSTER annottaion transfer to base object ###############################################

################################### RE-Subsetting Macrophage #############################################

# check cell type of macrophage
unique(macrophage$cluster_celltype)

# check ceöö type of base object 
unique(melanoma_wo_ud_macro_CD8$celltype_sc3)

# rename of celltype of macrophage by renaming undefined to macrophages
# now create a broader cell type category
cluster_celltype <- recode(macrophage$cluster_celltype,
                           "Monocytes" = "Monocytes",
                           "Myeloid" = "Myeloid",
                           "LATAM" = "LATAM",
                           "TAM" = "TAM",
                           "Antigen-Presenting Macrophages" = "Antigen-Presenting Macrophages",
                           "APOE+ cell" = "APOE+ macrophage",
                           "APC" = "APC",
                           "Undefined" = "Macrophages",
                           "IFNTAM" = "IFNTAM",
                           "Macrophages" ="Macrophages",
                           "Dendritic cells" = "Dendritic cells"
)

macrophage$cluster_celltype <- cluster_celltype

# Step 0: Create a mapping from `CD8`
annotation_mapping <- data.frame(
  cell_id = colnames(macrophage),                # Cell IDs from CD8
  cluster_celltype = macrophage$cluster_celltype # Annotations from CD8
)

# Step 1: Create a new column `celltype_sc3` and convert it to character
melanoma_wo_ud_macro_CD8$celltype_sc4 <- as.character(melanoma_wo_ud_macro_CD8$celltype_sc3)

# Step 2: Find matched indices
matched_indices <- match(colnames(melanoma_wo_ud_macro_CD8), annotation_mapping$cell_id)

# Step 3: Update `celltype_sc3` with the new annotations
melanoma_wo_ud_macro_CD8$celltype_sc4[!is.na(matched_indices)] <- as.character(
  annotation_mapping$cluster_celltype[matched_indices[!is.na(matched_indices)]]
)

# Step 4: (Optional) Convert `celltype_sc3` back to a factor
melanoma_wo_ud_macro_CD8$celltype_sc4 <- factor(melanoma_wo_ud_macro_CD8$celltype_sc4)

# Verify the updated annotations
table(melanoma_wo_ud_macro_CD8$celltype_sc4)

########################################################### Spatial graph analysis ###########################################################################################
# load filtered object

melanoma_wo_ud_macro_CD8 <- readRDS('/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/NEW_samples_NOV_2024/IMC/IRA_K15/melanoma_wo_ud_macro_CD8.rds')

#load object, images and masks
IMC <- readRDS('/Users/../IMC.rds')
images <- readRDS('/Users/../images.rds')
masks <- readRDS('/Users/../masks.rds')

# transfer cell annotation to IMC
melanoma_wo_ud_macro_CD8 <- buildSpatialGraph(melanoma_wo_ud_macro_CD8, img_id = "sample_id", type = "knn", k = 20)
melanoma_wo_ud_macro_CD8 <- buildSpatialGraph(melanoma_wo_ud_macro_CD8, img_id = "sample_id", type = "expansion", threshold = 20)
melanoma_wo_ud_macro_CD8 <- buildSpatialGraph(melanoma_wo_ud_macro_CD8, img_id = "sample_id", type = "delaunay", max_dist = 20)
colPairNames(melanoma_wo_ud_macro_CD8)

# define color of all cell types
# Define a color scheme for cell types
celltype_sc4_colors <- c(
  "Antigen-Presenting Macrophages" = "#F5F5DC", # Beige
  "APC" = "#F0F8FF", # Alice Blue
  "APOE+ macrophage" = "#FAEBD7", # Antique White
  "B cells" = "#FFFACD", # Lemon Chiffon
  "CAF" = "#E0FFFF", # Light Cyan
  "CD8 T cells" = "#00FF00", # Bright Green
  "Dendritic cells" = "#F5FFFA", # Mint Cream
  "EMT cells" = "#FDF5E6", # Old Lace
  "Endothelial progenitor" = "#F0FFF0", # Honeydew
  "Gymb+ CD8" = "#FF1493", # Bright Pink
  "Granulocytes" = "#F5DEB3", # Wheat
  "IFNTAM" = "#FFA500", # Bright Orange
  "Keratinocytes" = "#FFF5EE", # Seashell
  "LATAM" = "#FF0000", # Bright Red
  "Macrophages" = "#FFF8DC", # Cornsilk
  "Melanoma" = "#D3D3D3", # Light Gray
  "Monocytes" = "#FAFAD2", # Light Goldenrod Yellow
  "Myeloid" = "#FFE4E1", # Misty Rose
  "Perivascular CAF" = "#FFF0F5", # Lavender Blush
  "Stroma" = "#E6E6FA", # Light Lavender
  "T cells" = "#F0E68C", # Khaki
  "TIM+ CD8" = "#0000FF" # Bright Blue
)

# Define a color scheme for broad cell types
celltype_sc5_colors <- c(
  "Stroma" = "#ffb6c1", # Beige
  "B cells" = "#B2D8B2", # Alice Blue
  "IFNTAM" = "purple", # Antique White
  "Granulocytes" = "#FFFACD", # Lemon Chiffon
  "APC" = "#E0FFFF", # Light Cyan
  "CD8 T cells" = "#4e91fd", # Bright Green
  "Macrophages" = "#ff7b7b", # Mint Cream
  "Gymb+ CD8" = "#2c2cff", # Old Lace
  "LATAM" = "#ff0000", # Honeydew
  "EMT cells" = "#F5DEB3", # Wheat
  "Keratinocytes" = "#FFF5EE", # Seashell
  "Dendritic cells" = "#D3D3D3", # Light Gray
  "TIM+ CD8" = "#007300", # Light Goldenrod Yellow
  "Melanoma" = "#F0E68C", # Misty Rose
  "T cells" = "#40e0d0", # Khaki
  "Monocytes" = "#997a8d",
  "Endothelial progenitor" = "#b3ecec"
)


# Assign the color scheme to the metadata
metadata(melanoma_wo_ud_macro_CD8)$color_vectors$celltype_sc5 <- celltype_sc5_colors


# steinbock interaction graph 
plotSpatial(melanoma_wo_ud_macro_CD8[,melanoma_wo_ud_macro_CD8$sample_id == "25296303_5_5"], 
            node_color_by = "celltype_sc5", 
            img_id = "sample_id", 
            draw_edges = TRUE, 
            colPairName = "neighborhood", 
            nodes_first = FALSE,
            node_size_fix = 0.5,
            edge_color_fix = "white") + 
  scale_color_manual(values = metadata(melanoma_wo_ud_macro_CD8)$color_vectors$celltype_sc5) +
  ggtitle("steinbock interaction graph")

# inndividual expression
plotSpatial(melanoma_wo_ud_macro_CD8[,melanoma_wo_ud_macro_CD8$sample_id == "25296303_5_5"], 
            node_color_by = "GranzymeB", 
            assay_type = "exprs",
            img_id = "sample_id", 
            draw_edges = TRUE, 
            colPairName = "expansion_interaction_graph", 
            nodes_first = FALSE, 
            node_size_by = "area", 
            directed = FALSE,
            edge_color_fix = "grey") + 
  scale_size_continuous(range = c(0.1, 2)) +
  ggtitle("ZYMB expression")

# community
set.seed(230621)
melanoma_wo_ud_macro_CD8 <- detectCommunity(melanoma_wo_ud_macro_CD8, 
                                            colPairName = "neighborhood", 
                                            size_threshold = 10)

plotSpatial(melanoma_wo_ud_macro_CD8, 
            node_color_by = "spatial_community", 
            img_id = "sample_id", 
            node_size_fix = 0.5) +
  theme(legend.position = "none") +
  ggtitle("Spatial tumor communities") +
  scale_color_manual(values = rev(colors()))

# tumor- stroma communicty
# Define the mapping of cell types to compartments
architecture_mapping <- list(
  Tumor = c("Melanoma", "CAF", "Perivascular CAF", "EMT cells"),
  Stroma = c("Stroma", "Endothelial progenitor", "Keratinocytes", "Granulocytes"),
  Immune = c("Dendritic cells", "T cells", "TIM+ CD8", "Monocytes", "Myeloid",
             "CD8 T cells", "Gymb+ CD8", "TAM", "Antigen-Presenting Macrophages",
             "APOE+ macrophage", "APC", "B cells", "IFNTAM")
)

# Assign compartments based on cell type
melanoma_wo_ud_macro_CD8$architecture <- sapply(
  melanoma_wo_ud_macro_CD8$celltype_sc4,
  function(celltype) {
    if (celltype %in% architecture_mapping$Tumor) {
      return("Tumor")
    } else if (celltype %in% architecture_mapping$Stroma) {
      return("Stroma")
    } else if (celltype %in% architecture_mapping$Immune) {
      return("Immune")
    } else {
      return(NA)  # Handle any unexpected cell types
    }
  }
)

############################################### Interaction analysis ########################################

library(scales)
out <- testInteractions(melanoma_wo_ud_macro_CD8, 
                        group_by = "sample_id",
                        label = "celltype_sc5", 
                        colPairName = "neighborhood",
                        BPPARAM = SerialParam(RNGseed = 221029))

head(out)

# create heatmap
out %>% as_tibble() %>%
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE)) %>%
  ggplot() +
  geom_tile(aes(from_label, to_label, fill = sum_sigval)) +
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# now add a new column to minimize cell types by merging related cell types
# Create a lookup table with the original and new names

celltype_sc5 <- recode(melanoma_wo_ud_macro_CD8$celltype_sc4,
                       "Dendritic cells" = "Dendritic cells",
                       "T cells" = "T cells",
                       "TIM+ CD8" = "TIM+ CD8",
                       "Monocytes" = "Monocytes",
                       "Myeloid" = "Macrophages",
                       "CD8 T cells" = "CD8 T cells",
                       "Stroma" = "Stroma",
                       "Endothelial progenitor" = "Endothelial progenitor",
                       "Gymb+ CD8" = "Gymb+ CD8",
                       "CAF" = "Stroma",
                       "LATAM" = "LATAM",
                       "TAM" = "Macrophages",
                       "Melanoma" = "Melanoma",
                       "Antigen-Presenting Macrophages" = "Macrophages",
                       "APOE+ macrophage" = "Macrophages",
                       "Keratinocytes" = "Keratinocytes",
                       "Granulocytes" = "Granulocytes",
                       "APC" = "APC",
                       "Perivascular CAF" = "Stroma",
                       "Macrophages" = "Macrophages",
                       "B cells" = "B cells",
                       "EMT cells" = "EMT cells",
                       "IFNTAM" = "IFNTAM"
)

# Apply the mapping to the `celltype_sc4` column
melanoma_wo_ud_macro_CD8$celltype_sc5 <- celltype_sc5

# Verify the updated cell types
unique(melanoma_wo_ud_macro_CD8$celltype_sc5)
table(melanoma_wo_ud_macro_CD8$celltype_sc5)

######## subset only imp cell types ###########################
# Define the cell types to retain
celltypes_to_keep <- c("Dendritic cells", "TIM+ CD8", "Stroma", "Endothelial progenitor",
                       "Gymb+ CD8", "IFNTAM", "LATAM", "Melanoma", "B cells")

# Subset CD8_subset to keep only the specified cell types
filtered <- melanoma_wo_ud_macro_CD8[, melanoma_wo_ud_macro_CD8$celltype_sc5 %in% celltypes_to_keep]

# Verify the retained cell types
unique(filtered$celltype_sc5)


######################## Subsetting based on response #######################################

# subset based on response column 
non_responder <- melanoma_wo_ud_macro_CD8[, melanoma_wo_ud_macro_CD8$Response == "Non-Responder"] 
responder <- melanoma_wo_ud_macro_CD8[, melanoma_wo_ud_macro_CD8$Response == "Responder"]

# subset based on response column 
non_responder <- filtered[, filtered$Response == "Non-Responder"] 
responder <- filtered[, filtered$Response == "Responder"]


#### Responder
out <- testInteractions(responder, 
                        group_by = "sample_id",
                        label = "celltype_sc5", 
                        colPairName = "expansion_interaction_graph",
                        BPPARAM = SerialParam(RNGseed = 221029))

head(out)

# create heatmap
out %>% as_tibble() %>%
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE)) %>%
  ggplot() +
  geom_tile(aes(from_label, to_label, fill = sum_sigval)) +
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#### non-Responder
out <- testInteractions(non_responder, 
                        group_by = "sample_id",
                        label = "celltype_sc5", 
                        colPairName = "neighborhood",
                        BPPARAM = SerialParam(RNGseed = 221029))

head(out)

# create heatmap
out %>% as_tibble() %>%
  group_by(from_label, to_label) %>%
  summarize(sum_sigval = sum(sigval, na.rm = TRUE)) %>%
  ggplot() +
  geom_tile(aes(from_label, to_label, fill = sum_sigval)) +
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##################################### Final Image Analysis ###############################################
#load object, images and masks
IMC <- readRDS('/Users/../IMC.rds')
images <- readRDS('/Users/../images.rds')
masks <- readRDS('/Users/../masks.rds')



# Sample images
set.seed(220517)
# Predefined IDs you want to include

# non_res_id
non_res_id <- "20_2_2"

# Remaining IDs (excluding the specific ones)
res_id <- "03_5_5"

# For Non Responder Patients
#cur_id <- sample(unique(IMC$sample_id), 9)
cur_images <- images[names(images) %in% non_res_id] # non_responder patients
cur_masks <- masks[names(masks) %in% non_res_id] # non_responder patients

# For Responder Patients
cur_images <- images[names(images) %in% res_id] # responder patients
cur_masks <- masks[names(masks) %in% res_id] # responder patients

#Pixel visualization APOE+ TERM2
plotPixels(cur_images, 
           colour_by = c("APOE","TREM2", "CD8a"),
           bcg = list(APOE = c(0, 2, 1),
                      TREM2 = c(0, 10, 1),
                      CD8a = c(0, 5, 1)))

#Pixel visualization APOE+ TERM2
plotPixels(cur_images, 
           colour_by = c("GranzymeB","CCL4", "CD8a"),
           bcg = list(GranzymeB = c(0, 5, 1),
                      CCL4 = c(0, 1, 1),
                      CD8a = c(0, 5, 1)))
