library(Seurat)
library(ggplot2)

all <- readRDS('Source Data II_Seurat Files/all.rds')

DimPlot(all)

my <- readRDS('my.rds')

DefaultAssay(my) <- "RNA"
DimPlot(my, reduction = "umap")

# Subset to exclude subcutaneous ("sc") tissue
adipose_filtered <- subset(my, subset = tissue != "sc")

DimPlot(adipose_filtered, reduction = "umap")

# Create named vector of annotations
cluster_labels <- c(
  "0"  = "M2-like macrophage (MRC1⁺, F13A1⁺)",
  "1"  = "M2-like macrophage (MRC1⁺, RBPJ⁺)",
  "2"  = "LAM (TREM2⁺, CD9⁺, FABP4⁺)",
  "3"  = "cDC2 (CD1c⁺, CLEC10A⁺, HLA-DR⁺)",
  "4"  = "LYVE1⁺ M2-like macrophage (tissue support)",
  "5"  = "Non-classical monocyte (FCGR3A⁺, TCF7L2⁺)",
  "6"  = "M1-/M2-like macrophage (CD11c⁺, CTSB⁺)",
  "7"  = "Vasculature-associated macrophage (FN1⁺, MARCO⁺)",
  "8"  = "M2-like macrophage (ECM remodeling, TIMP1⁺)",
  "9"  = "M2-like macrophage (TGFBI⁺, ABCA1⁺, PPARG⁺)",
  "10" = "Metabolic macrophage (MSR1⁺, CD36⁺, PLIN2⁺)",
  "11" = "C1Qhi M2-like macrophage (SEPP1⁺, RNASE1⁺)",
  "12" = "M2-like macrophage (C3⁺, STAT3⁺, FKBP5⁺)",
  "13" = "Autophagy-high M2 macrophage (MERTK⁺, SQSTM1⁺)",
  "14" = "Classical monocyte (S100A8⁺, CD14⁺)",
  "15" = "MOX macrophage (HMOX1⁺, IL1B⁺, NRF2⁺)"
)

# Convert Idents to character so they match the names in cluster_labels
cluster_ids <- as.character(Idents(adipose_filtered))

# Map labels to each cell
myeloid_labels <- cluster_labels[cluster_ids]

# Assign names to match Seurat object barcodes
names(myeloid_labels) <- colnames(adipose_filtered)

# Add as metadata column
adipose_filtered$myeloid_cluster_label <- myeloid_labels

DimPlot(adipose_filtered, group.by = "myeloid_cluster_label", label = TRUE, repel = TRUE) + NoLegend()

# Set identities to clusters if not already
Idents(adipose_filtered) <- "seurat_clusters"  # or keep as current if already correct

adipose_filtered <- ScaleData(adipose_filtered, verbose = FALSE)

# Violin plot for all four genes
VlnPlot(
  object = adipose_filtered,
  features = c("TREM2", "CD9", "FABP4", "CD68"),
  slot = "scale.data",
  pt.size = 0,                  # hides individual dots (optional)
  stack = TRUE,                 # stack vertically
  flip = TRUE                   # clusters on y-axis
)



macro <- readRDS('feldman_seurat_UCell_SC.rds')

################################## Integration of macro and adipose_filtered ###########################
# normalize data for adipose filter
adipose_filtered <- NormalizeData(
  object = adipose_filtered,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000
)


#step 1
# rename seurat cluster accroding to dataset names
#check colnames
colnameLAM <- colnames(adipose_filtered@meta.data)
adipose_filtered$LAM.cluster <- as.character(Idents(adipose_filtered))
colnameTAM <- colnames(macro@meta.data)
# Copy seurat_clusters to new column
macro$TAM.cluster <- macro$seurat_clusters

# Optional: remove the original column (if desired)
macro@meta.data$seurat_clusters <- NULL

#step 2 
# add a column with dataset name
# Add a new column to meta.data in seurat_myeloid
adipose_filtered@meta.data$dataset <- "LAM"

# Add a new column to meta.data in seurat_monocyte
macro@meta.data$dataset <- "TAM"

# Verify the changes
head(adipose_filtered@meta.data)
head(macro@meta.data)

# create a list of seurat object 
seurat_list <- list(adipose_filtered, macro)

# Name the elements of the list
names(seurat_list) <- c("LAM", "TAM")

######################################################################################################
########################## Seurat Anchoring ##########################################################
######################################################################################################

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = seurat_list)
seurat_list <- lapply(X = seurat_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})


# anchor based integration
LAM_TAM_anchors <- FindIntegrationAnchors(object.list = seurat_list, reference = 2, anchor.features = features, reduction = "rpca")

# this command creates an 'integrated' data assay
LAM_TAM_combined <- IntegrateData(anchorset = LAM_TAM_anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(LAM_TAM_combined) <- "integrated"


# Run the standard workflow for visualization and clustering
LAM_TAM_combined <- ScaleData(LAM_TAM_combined, verbose = FALSE)
LAM_TAM_combined <- RunPCA(LAM_TAM_combined, npcs = 20, verbose = FALSE)
LAM_TAM_combined <- RunUMAP(LAM_TAM_combined, reduction = "pca", dims = 1:20)
LAM_TAM_combined <- FindNeighbors(LAM_TAM_combined, reduction = "pca", dims = 1:20)
LAM_TAM_combined <- FindClusters(LAM_TAM_combined, resolution = 0.6)

# Visualization
p1 <- DimPlot(LAM_TAM_combined, reduction = "umap", group.by = "dataset")

p1

################ highlight cluster 2 of TAM and LAM ######################
# Create a new metadata column for highlighting
LAM_TAM_combined$highlight <- "Other"  # default label for all cells

# Assign labels based on cluster membership
LAM_TAM_combined$highlight[
  LAM_TAM_combined$TAM.cluster == 2
] <- "TAM cluster 2"

LAM_TAM_combined$highlight[
  LAM_TAM_combined$LAM.cluster == 2
] <- "LAM cluster 2"

# Plot with custom colors
DimPlot(
  LAM_TAM_combined,
  reduction = "umap",
  group.by = "highlight",
  cols = c(
    "Other" = "gray80",
    "TAM cluster 2" = "#1f78b4",  # blue
    "LAM cluster 2" = "#e31a1c"   # red
  ),
  pt.size = 1
) + ggtitle("UMAP: Highlighting Cluster 2 in TAM and LAM")


FeaturePlot(
  object = LAM_TAM_combined,
  features = c("APOE", "TREM2", "CD9", "CD36"),
  reduction = "umap",
  cols = c("lightgray", "red"),  # adjust colors as needed
  pt.size = 1,
  ncol = 2,
  order = TRUE
)

# Use actual expression data
DefaultAssay(LAM_TAM_combined) <- "RNA"  # or "SCT"

FeaturePlot(
  LAM_TAM_combined,
  features = c("APOE", "TREM2", "CD9", "CD36"),
  reduction = "umap",
  slot = "data",           # default, log-normalized expression
  cols = c("lightgray", "red"),
  pt.size = 1,
  ncol = 2,
  order = TRUE
)



#################### Combined cluster #############################
# Initialize with "NA"
LAM_TAM_combined$combined.cluster <- NA

# Assign labels from TAM.cluster (e.g., 1 → TAM1)
tam_mask <- !is.na(LAM_TAM_combined$TAM.cluster)
LAM_TAM_combined$combined.cluster[tam_mask] <- paste0("TAM", LAM_TAM_combined$TAM.cluster[tam_mask])

# Assign labels from LAM.cluster (e.g., 1 → LAM1)
lam_mask <- !is.na(LAM_TAM_combined$LAM.cluster)
LAM_TAM_combined$combined.cluster[lam_mask] <- paste0("LAM", LAM_TAM_combined$LAM.cluster[lam_mask])

table(LAM_TAM_combined$combined.cluster)

DimPlot(
  LAM_TAM_combined,
  reduction = "umap",
  group.by = "combined.cluster",
  label = TRUE,
  repel = TRUE,
  pt.size = 1
) + ggtitle("UMAP: Combined LAM and TAM Clusters")

LAM_TAM_combined <- NormalizeData(
  object = LAM_TAM_combined,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

DefaultAssay(LAM_TAM_combined) <- "RNA"

FeaturePlot(
  LAM_TAM_combined,
  features = c("APOE", "TREM2", "CD9", "CD36"),
  reduction = "umap",
  slot = "data",
  cols = c("lightgray", "red"),
  pt.size = 1,
  ncol = 2,
  order = TRUE
)

############################ Violin plot #####################################


VlnPlot(
  object = LAM_TAM_combined,
  features = c("APOE", "TREM2", "CD36", "CD9"), 
  group.by = "combined.cluster",
  slot = "data",             # use log-normalized expression
  pt.size = 0,               # hides single-cell points (optional)
  flip = TRUE                # horizontal violins
)








######################## End Violin #########################################

# save LAM_TAM_combined
saveRDS(LAM_TAM_combined, 'LAM_TAM_combined.RDS')

# load
# save LAM_TAM_combined
LAM_TAM_combined <- readRDS('LAM_TAM_combined.RDS')

VlnPlot(
  object = LAM_TAM_combined,
  features = c("CD163", "CD63", "CD68", "MMP9", "MMP7"), 
  group.by = "combined.cluster",
  slot = "data",             # use log-normalized expression
  pt.size = 0,               # hides single-cell points (optional)
  flip = TRUE                # horizontal violins
)


########################## DEfine Custom Colors ###################
# Define LAM colors (cool tones)
lam_colors <- setNames(
  c(
    "#1f78b4", "#33a02c", "#6a3d9a", "#a6cee3", "#b2df8a",
    "#cab2d6", "#8dd3c7", "#80b1d3", "#bebada", "#b3de69",
    "#fccde5", "#bc80bd", "#ccebc5", "#ffed6f", "#66c2a5",
    "#a6d854"
  ),
  paste0("LAM", 0:15)
)

# Define TAM colors (warm tones)
tam_colors <- setNames(
  c(
    "#e31a1c", "#fb9a99", "#ff7f00", "#fdbf6f", "#b15928",
    "#fdae61", "#d95f02", "#e78ac3", "#d73027"
  ),
  paste0("TAM", 0:8)
)

# Combine into one vector
combined_colors <- c(tam_colors, lam_colors)

# UMAP plot
DimPlot(
  LAM_TAM_combined,
  reduction = "umap",
  group.by = "combined.cluster",
  label = TRUE,
  repel = TRUE,
  pt.size = 1,
  cols = combined_colors
) + ggtitle("UMAP: Combined LAM and TAM Clusters")


########################## transcriptional concordance  #############################################
library(Seurat)
library(pheatmap)

# 1. Use log-normalized expression (slot = "data")
avg_expr <- AverageExpression(
  LAM_TAM_combined,
  group.by = "combined.cluster",
  assays = "RNA",
  slot = "data"
)[["RNA"]]  # genes x clusters

# 2. Rank gene expression within each cluster
# Each column gets ranked independently (higher expression → higher rank)
rank_expr <- apply(avg_expr, 2, function(x) rank(x, ties.method = "average"))

# 3. Transpose to get clusters as rows
cluster_rank_expr <- t(rank_expr)  # [clusters x genes]

# 4. Define TAM and LAM cluster names
tam_names <- grep("^TAM", rownames(cluster_rank_expr), value = TRUE)
lam_names <- grep("^LAM", rownames(cluster_rank_expr), value = TRUE)

cor_rank <- cor(
  t(cluster_rank_expr[tam_names, ]),  # genes x TAM clusters
  t(cluster_rank_expr[lam_names, ]),  # genes x LAM clusters
  method = "pearson"
)

# 6. Visualize heatmap
pheatmap(
  mat = cor_rank,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",
  main = "TAM vs LAM Cluster Correlation (Gene Rank)",
  fontsize = 10,
  color = colorRampPalette(c("blue", "white", "red"))(100)
)


####################################################################################################################
####################################################################################################################
################################## DECOUPLER single cells ##########################################################
####################################################################################################################
####################################################################################################################
## We load the required packages
library(Seurat)
library(decoupleR)
library(openxlsx)
# Only needed for data handling and plotting
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)

############# load net (gene set) data #######################

net <- read.xlsx('SPATIAL_METABOLISM/NET.xlsx')

# Extract the normalized log-transformed counts
mat <- as.matrix(LAM_TAM_combined@assays$RNA@data)

# Run mlm
acts <- decoupleR::run_ulm(mat = mat, 
                           net = net, 
                           .source = 'source', 
                           .target = 'target',
                           .mor = 'weight', 
                           minsize = 5)
acts

# Extract mlm and store it in pathwaysmlm in data
LAM_TAM_combined[['pathwaysulm']] <- acts %>%
  tidyr::pivot_wider(id_cols = 'source', 
                     names_from = 'condition',
                     values_from = 'score') %>%
  tibble::column_to_rownames(var = 'source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
Seurat::DefaultAssay(object = LAM_TAM_combined) <- "pathwaysulm"

# Scale the data
LAM_TAM_combined <- Seurat::ScaleData(LAM_TAM_combined)

LAM_TAM_combined@assays$pathwaysulm@data <- LAM_TAM_combined@assays$pathwaysulm@scale.data

p1 <- Seurat::DimPlot(LAM_TAM_combined, 
                      reduction = "umap", 
                      label = TRUE, 
                      pt.size = 0.5) + 
  Seurat::NoLegend() + 
  ggplot2::ggtitle('Cell types')

colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])

p2 <- Seurat::FeaturePlot(LAM_TAM_combined, features = c("Lipid-Uptake-and-Transport")) + 
  ggplot2::scale_colour_gradient2(low = colors[1], mid = 'white', high = colors[2]) +
  ggplot2::ggtitle('Lipid-Uptake-and-Transport')

p <- p1 | p2
p

#Define features to plot
features_to_plot <- c("Cholesterol-Metabolism", 
                      "Fatty-Acid-Oxidation", 
                      "Fatty-Acid-Synthesis",
                      "Glycolysis",
                      "Lipid-Uptake-and-Transport",
                      "Lipolysis",
                      "Oxidative-Phosphorylation",
                      "ROS-Production")  # update with your actual feature names

# 3. Define custom diverging colors
colors <- c("blue", "red")  # low = blue, high = red

feature_plots <- lapply(features_to_plot, function(f) {
  FeaturePlot(LAM_TAM_combined, 
              features = f, 
              reduction = "umap", 
              pt.size = 0.2,       # smaller points
              order = TRUE         # plot high values on top
  ) +
    scale_colour_gradient2(
      low = colors[1], 
      mid = "white", 
      high = colors[2], 
      midpoint = 0
    ) +
    ggtitle(f)
})

# Combine and display
wrap_plots(feature_plots, ncol = 3)

########################## HEATMAP ##################################

data <- LAM_TAM_combined

Idents(data) <- data$combined.cluster

# Extract activities from object as a long dataframe
df <- t(as.matrix(data@assays$pathwaysulm@data)) %>%
  as.data.frame() %>%
  dplyr::mutate(cluster = Seurat::Idents(data)) %>%
  tidyr::pivot_longer(cols = -cluster, 
                      names_to = "source", 
                      values_to = "score") %>%
  dplyr::group_by(cluster, source) %>%
  dplyr::summarise(mean = mean(score))

# Transform to wide matrix
top_acts_mat <- df %>%
  tidyr::pivot_wider(id_cols = 'cluster', 
                     names_from = 'source',
                     values_from = 'mean') %>%
  tibble::column_to_rownames(var = 'cluster') %>%
  as.matrix()

top_acts_df <- as.data.frame(top_acts_mat)
write.xlsx(top_acts_df, 'top_acts.xlsx', rowNames = TRUE)

top_acts_df <- read.xlsx('top_acts.xlsx', rowNames = TRUE)

top_acts_mat <- as.matrix(top_acts_df)

# Color scale
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)

my_breaks <- c(seq(-1.25, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.05, 1.25, length.out = floor(100 / 2)))

# Plot
pheatmap::pheatmap(mat = top_acts_mat,
                   color = colors.use,
                   border_color = "white",
                   breaks = my_breaks,
                   cellwidth = 20,
                   cellheight = 20,
                   treeheight_row = 20,
                   treeheight_col = 20) 


############################# Lipid genes dot plot ##################################
library(viridis)
library(RColorBrewer)

# laod gene list
goi <- read.xlsx('Gene_list_UKF_sample.xlsx', sheet = 2)

library(tidyr)
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)


Idents(LAM_TAM_combined) <- LAM_TAM_combined$combined.cluster

# 1. Convert goi to long format and remove NA
goi_long <- pivot_longer(goi, cols = everything(), names_to = "Pathway", values_to = "Gene") %>%
  filter(!is.na(Gene))

DefaultAssay(LAM_TAM_combined) <- "RNA"

# Check available genes
available_genes <- rownames(LAM_TAM_combined)

# Filter for valid gene names only
goi_long_valid <- goi_long %>%
  filter(Gene %in% available_genes)

# 2. Create named list of features grouped by pathway
goi_list <- goi_long_valid %>%
  group_by(Pathway) %>%
  summarise(genes = list(Gene)) %>%
  deframe()  # Named list: Pathway -> vector of genes

all_genes <- unique(unlist(goi_list))

DotPlot(LAM_TAM_combined, features = all_genes) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Gene DotPlot grouped by Pathway")

########################### PRE AND POST Enrichment comparison #########################################
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# 1. Set the default assay to pathwaysulm
DefaultAssay(macro) <- "pathwaysulm"

# 2. Get pathway score names (assumed stored in assay rows)
pathway_scores <- rownames(macro)  # e.g., Lipid-Uptake-and-Transport, etc.

#rename
macro@meta.data <- macro@meta.data %>%
  rename(
    Treatment_Status = `characteristics..patinet.ID..Pre.baseline..Post..on.treatment.`,
    Response = `characteristics..response`
  )

# 3. Fetch scores + metadata
df <- FetchData(macro, vars = c(pathway_scores, "seurat_clusters", "Treatment_Status", "Response"))

df <- df %>%
  filter(seurat_clusters == 2)

df$Treatment_Status <- sub("_.*", "", df$Treatment_Status)

# Define selected pathways
selected_pathways <- c("Cholesterol-Metabolism", 
                       "Lipid-Uptake-and-Transport", 
                       "Lipolysis")

# convert to long
# Reshape to long format
df_long <- df %>%
  pivot_longer(cols = all_of(selected_pathways), 
               names_to = "Pathway", 
               values_to = "Score")


df_long$Treatment_Status <- factor(df_long$Treatment_Status, levels = c("Pre", "Post"))

#plot
# Violin plot with p-values
ggplot(df_long, aes(x = Treatment_Status, y = Score, fill = Treatment_Status)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.6) +
  facet_grid(rows = vars(Response), cols = vars(Pathway), scales = "free_y") +
  scale_fill_manual(values = c("Pre" = "#377eb8", "Post" = "#e41a1c")) +
  theme_minimal(base_size = 13) +
  labs(title = "Pre vs Post Pathway Scores in Cluster 2",
       x = "Treatment Status", 
       y = "Pathway Score") +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("Pre", "Post")),
    label = "p.adj",
    method.args = list(p.adjust.method = "BH"),
    label.y = 1.1  # optional: raise label position
  )

############################## Responder and non-Responder in Pre and Post ############################
library(ggplot2)
library(ggpubr)

# 1. Ensure proper factor levels for consistent ordering
df_long$Response <- factor(df_long$Response, levels = c("Non-responder", "Responder"))
df_long$Treatment_Status <- factor(df_long$Treatment_Status, levels = c("Pre", "Post"))

# 2. Plot Responder vs Non-responder comparison, split by Pre/Post
ggplot(df_long, aes(x = Response, y = Score, fill = Response)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.6) +
  facet_grid(rows = vars(Treatment_Status), cols = vars(Pathway), scales = "free_y") +
  scale_fill_manual(values = c("Responder" = "#377eb8", "Non-responder" = "#e41a1c")) +
  theme_minimal(base_size = 13) +
  labs(title = "Responder vs Non-responder by Pathway and Treatment",
       x = "Response Group", 
       y = "Pathway Score") +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("Non-responder", "Responder")),
    label = "p.adj",
    method.args = list(p.adjust.method = "BH"),
    label.y = 1.1  # adjust based on data range
  )

################## ONLY res and nonRes ###############################
# Plot across responders vs non-responders, for each pathway
ggplot(df_long, aes(x = Response, y = Score, fill = Response)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.6) +
  facet_wrap(~ Pathway, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("Responder" = "#377eb8", "Non-responder" = "#e41a1c")) +
  theme_minimal(base_size = 13) +
  labs(title = "Responder vs Non-responder (All Treatment Statuses)",
       x = "Response Group", 
       y = "Pathway Score") +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("Non-responder", "Responder")),
    label = "p.adj",
    method.args = list(p.adjust.method = "BH"),
    label.y = 1.1  # adjust this based on score range
  )

############################### LAM vs. LA-TAM DEG #######################################
all_genes <- rownames(LAM_TAM_combined)

# Expanded patterns for filtering
remove_patterns <- c(
  "^MT-", "^RPL", "^RPS",              # mitochondrial & ribosomal
  "^MIR", "^MIRLET",                   # microRNAs
  "^RNU", "^SNORD", "^SCARNA",         # small RNAs
  "^HIST",                             # histone genes
  "^RP[0-9]+-", "^CTD-", "^AC", "^AL", # unannotated loci
  "^LOC", "^AP", "^FAM",               # orphan/family genes
  "MALAT1", "NEAT1", "RNASEK-C17orf49",# known lncRNAs
  "orf", "P[0-9]+$",                   # ORFs and pseudogenes
  "C[0-9]+orf[0-9]+", "KIAA"           # C1orf85-type, KIAA genes
)

# Combine regex
remove_regex <- paste(remove_patterns, collapse = "|")

# Identify and remove
genes_to_remove <- grep(remove_regex, all_genes, value = TRUE, ignore.case = TRUE)
genes_to_keep <- setdiff(all_genes, genes_to_remove)

# Subset object
LAM_TAM_combined <- subset(LAM_TAM_combined, features = genes_to_keep)

unique(Idents(LAM_TAM_combined))
Idents(LAM_TAM_combined) <- LAM_TAM_combined$combined.cluster

# DEG
LAM_vs_TAM <- FindMarkers(LAM_TAM_combined, ident.1 = "LAM2", ident.2 = "TAM2")

# view results
head(LAM_vs_TAM)

library(openxlsx)

write.xlsx(LAM_vs_TAM, 'LAM_TAM_DEG.xlsx', rowNames = TRUE)

##########################################Enrichment ####################################################
# Load libraries
library(dplyr)
library(STRINGdb)

#load
LAM_vs_TAM <- read.xlsx('LAM_TAM_DEG.xlsx')

# ---------------------------
# 1. Filter significant DEGs
# ---------------------------
deg_filtered <- LAM_vs_TAM %>%
  filter(p_val_adj < 0.01)

# ---------------------------
# 2. Split into UP and DOWN
# ---------------------------
top_up <- deg_filtered %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 200)

top_down <- deg_filtered %>%
  arrange(avg_log2FC) %>%
  slice_head(n = 200)


# 
write.xlsx(top_down, 'LAM_TopDOWN200.xlsx')


# ---------------------------
# 3. Prepare gene lists
# ---------------------------
genes_up <- rownames(top_up)
genes_down <- rownames(top_down)


# 4. Initialize STRINGdb (choose organism: 9606 for human)
string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400, input_directory = "")

# ---------------------------
# 5. Map and run enrichment
# ---------------------------

## UP
mapped_up <- string_db$map(data.frame(gene = genes_up), "gene", removeUnmappedRows = TRUE)
enrich_up <- string_db$get_enrichment(mapped_up$STRING_id)

## DOWN
mapped_down <- string_db$map(data.frame(gene = genes_down), "gene", removeUnmappedRows = TRUE)
enrich_down <- string_db$get_enrichment(mapped_down$STRING_id)

# ---------------------------
# 6. View results
# ---------------------------
head(enrich_up)
head(enrich_down)

write.xlsx(enrich_up, 'LAM_UP2.xlsx')
write.xlsx(enrich_down, 'TAM_UP2.xlsx')

#load final version
enrichment_df <- read.xlsx('TOP5_GO_KRGG_LAM_UP_DOWN.xlsx')

# Select top 20 per group (or change to top 10 if you prefer)
top_terms <- enrichment_df %>%
  group_by(Regulation) %>%
  slice_max(order_by = abs(signed_logFDR), n = 20, with_ties = FALSE) %>%
  ungroup()

# Sort DOWN from most negative to least, and UP from most positive to least
top_terms <- top_terms %>%
  arrange(Regulation, desc(signed_logFDR)) %>%
  mutate(term_full = factor(term_full, levels = term_full))

library(ggplot2)

ggplot(top_terms, aes(x = term_full, y = signed_logFDR, fill = Regulation)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("UP" = "#2c7bb6", "DOWN" = "#d7191c")) +
  labs(
    x = NULL,
    y = expression(paste("Signed ", -log[10], "(FDR)")),
    title = "LAM_vs_LA-TAM",
    fill = "Regulation"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )
