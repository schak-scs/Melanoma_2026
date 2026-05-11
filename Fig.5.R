# load dblist
dbList <- "/Users/sajibchakraborty/Documents/proteomics_genesets/dbList_final2.rds"
dbList <- readRDS("/Users/sajibchakraborty/Documents/proteomics_genesets/dbList_final2.rds")
names(dbList)


############################# EXtract terms associated with METABOLISM ########################
# Filter KEGG terms containing "metabolism"
kegg_terms <- names(dbList$KEGG)
kegg_metabolism_idx <- grepl("metabolism", kegg_terms, ignore.case = TRUE)
kegg_metabolism_list <- dbList$KEGG[kegg_metabolism_idx]

# Convert to data frame
kegg_metabolism_df <- data.frame(
  Term = names(kegg_metabolism_list),
  Genes = sapply(kegg_metabolism_list, function(genes) paste(genes, collapse = ", ")),
  stringsAsFactors = FALSE
)

# Filter REACTIOME terms containing "metabolism"
# For REACTIOME
reactome_terms <- names(dbList$REACTOME)
reactome_metabolism_idx <- grepl("metabolism", reactome_terms, ignore.case = TRUE)
reactome_metabolism_list <- dbList$REACTOME[reactome_metabolism_idx]

reactome_metabolism_df <- data.frame(
  Term = names(reactome_metabolism_list),
  Genes = sapply(reactome_metabolism_list, function(genes) paste(genes, collapse = ", ")),
  stringsAsFactors = FALSE
)

# merge KEGG and REACTOME data frames
# Add a source column to each to track origin
kegg_metabolism_df$Source <- "KEGG"
reactome_metabolism_df$Source <- "REACTOME"

# Merge the two data frames
metabolism_combined_df <- rbind(kegg_metabolism_df, reactome_metabolism_df)

# View the merged result
head(metabolism_combined_df)

write.xlsx(metabolism_combined_df, "/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATIAL_METABOLISM/KEGG_REACTOME.xlsx")

################################### LAOD SPATIAL data ##################################
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

#load data
UKF3 <- readRDS('/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/NEW_samples_NOV_2024/Spatial_data_Nov_2024/processed/S3_NonRes_NEW/S3_nonres_clustering.RDS')
UKF5 <- readRDS('/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/NEW_samples_NOV_2024/Spatial_data_Nov_2024/processed/S5_Res_OLDS3/S5_spata_object_clu_GSEA.RDS')
UKF6 <- readRDS('/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/NEW_samples_NOV_2024/Spatial_data_Nov_2024/processed/S6_Res_NEW/S6_res_clustering.RDS')
UKF1 <- readRDS('/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATIAL_CELL_TO_CELL_COMM/UKF1new.rds')
UKF2 <- readRDS('/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATIAL_CELL_TO_CELL_COMM/UKF2new.rds')
UKF4 <- readRDS('/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATIAL_CELL_TO_CELL_COMM/UKF4new.rds')

#rename
#object <- UKF1
#object <- UKF2
#object <- UKF3
#object <- UKF4
#object <- UKF5
#object <- UKF6

# normalize if needed
UKF1 <- normalizeCounts(UKF1, method = "CLR")
UKF2 <- normalizeCounts(UKF2, method = "CLR")
UKF3 <- normalizeCounts(UKF3, method = "CLR", overwrite = TRUE)
UKF4 <- normalizeCounts(UKF4, method = "CLR")
UKF5 <- normalizeCounts(UKF5, method = "CLR")
UKF6 <- normalizeCounts(UKF6, method = "CLR")


# normalize if need be
getMatrixNames(object)
object <- normalizeCounts(object, method = "CLR")
object <- activateMatrix(object, mtr_name = "counts")

# results are immediately stored in the objects feature data
getGroupingOptions(object)

# plot bayes space cluster
plotSurface(object, color_by = "Lloyd_k6", pt_clrp = "uc")

# run DEA analysis
# bayes_space
object <- runDEA(object = object, across = "bayes_space", method = "wilcox")

# extract the complete data.frame
dea_df <- 
  getDeaResultsDf(
    object = object, 
    across = "bayes_space"
  )

head(dea_df)

# extract only significant genes
# e.g. top 10 genes for histology area 'tumor' 
dea_df_filter <- getDeaResultsDf(
  object = object, 
  across = "bayes_space",
  max_adj_pval = 0.05 # pval must be lower or equal than 0.01
)

######################### Gene Set enrichment for metabolism #################################
# Create a named list of gene sets
metabolism_gene_sets <- split(
  strsplit(metabolism_combined_df$Genes, ",\\s*") |> 
    lapply(trimws),
  metabolism_combined_df$Term
)

library(dplyr)
library(purrr)

# Background gene universe
gene_universe <- unique(dea_df_filter$gene)

# Function to test enrichment for one cluster
run_enrichment <- function(cluster_id) {
  
  # Get genes in this cluster
  cluster_genes <- dea_df_filter %>%
    filter(bayes_space == cluster_id) %>%
    pull(gene) %>%
    unique()
  
  # Test each metabolism gene set
  enrichment_results <- map_dfr(names(metabolism_gene_sets), function(term) {
    term_genes <- metabolism_gene_sets[[term]]
    
    # Build contingency table
    overlap <- length(intersect(cluster_genes, term_genes))
    only_in_cluster <- length(setdiff(cluster_genes, term_genes))
    only_in_term <- length(setdiff(term_genes, cluster_genes))
    neither <- length(setdiff(gene_universe, union(cluster_genes, term_genes)))
    
    contingency <- matrix(c(overlap, only_in_term, only_in_cluster, neither), nrow = 2)
    
    # Fisher's test
    test <- fisher.test(contingency, alternative = "greater")
    
    tibble(
      Cluster = cluster_id,
      Term = term,
      Overlap = overlap,
      SetSize = length(term_genes),
      ClusterSize = length(cluster_genes),
      PValue = test$p.value
    )
  })
  
  enrichment_results
}

# Run for all clusters
cluster_ids <- unique(dea_df_filter$bayes_space)
all_enrichment <- map_dfr(cluster_ids, run_enrichment)

# Adjust p-values
all_enrichment <- all_enrichment %>%
  group_by(Cluster) %>%
  mutate(adj_pval = p.adjust(PValue, method = "BH")) %>%
  ungroup()

# view
all_enrichment %>%
  filter(adj_pval < 0.05) %>%
  arrange(Cluster, adj_pval)


############################# surface plot for individual genes ###################

S1P <- c("SPHK1", "SPHK2", "SGPP1", "SGPP2", "S1PR1", "S1PR2", "S1PR3", "S1PR4", "S1PR5", "SGPL1")

spata.genes <- getGenes(object)
GOI <- intersect(S1P, spata.genes)

# compare gene expression on the surface
plotSurfaceComparison(
  object = object, 
  color_by = GOI, 
  pt_clrsp = "Greens 3", 
  display_image = TRUE, 
  #smooth = TRUE, 
  alpha_by = TRUE
) 

# compare gene expression on the surface
plotSurfaceComparison(
  object = object, 
  color_by = GOI,
  pt_size = 2,
  smooth_span = 0.2,
  pt_clrsp = "Heat", 
  display_image = FALSE, 
  #smooth = TRUE, 
  alpha_by = TRUE
) 

# ssGSEA

# gene set enrichment
# GSEA
# DEFINE TAM SIGNATURES
geneMat_TAM <- read.xlsx(file.path('/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATIAL_METABOLISM/LAM_Meta.xlsx'))
geneList_TAM <- lapply(1:ncol(geneMat_TAM), function(i) as.character(geneMat_TAM[, i]))
names(geneList_TAM) <- colnames(geneMat_TAM)
geneList_TAM <- lapply(geneList_TAM, function(i) unique(i[!is.na(i)]))


# GENES
#spata.genes.meta <- getGeneMetaData(spata_obj)
spata.genes <- getGenes(object)
geneList_TAM <- lapply(geneList_TAM, intersect, y = spata.genes)


# Loop through gene sets and add to the SPATA2 object
for (gs_name in names(geneList_TAM)) {
  genes <- geneList_TAM[[gs_name]]
  
  object <- addGeneSet(
    object = object,
    genes = genes,
    name = gs_name,
    class = "TAM",        # use a consistent custom class label
    overwrite = TRUE,
    check = TRUE
  )
}

##############
# PLOT SURFACE

genesetS2 <- c(
               "TAM_Lipid_Uptake_and_Transport",
               "TAM_Fatty_Acid_Metabolism_General",
               "TAM_Fatty_Acid_Oxidation",
               "TAM_Fatty_Acid_Synthesis",
               "TAM_Lipolysis",
               "TAM_Cholesterol_Metabolism",
               "TAM_Oxidative_Phosphorylation",
               "TAM_Glycolysis",
               "TAM_ROS_Production")

# surface plot 

# open application to obtain a list of plots
#plots <- plotSurfaceInteractive(object = spata_obj)

# ssGSEA
Spata_GS <- plotSurfaceComparison(object = object, 
                                  color_by = genesetS2,
                                  method_gs = "plage",
                                  smooth = TRUE,
                                  pt_clrsp = "inferno",
                                  smooth_span = 0.2,
                                  pt_size = 2,
                                  display_image = FALSE
)

Spata_GS


# plot bayes space cluster
plotSurface(object, color_by = "bayes_space", pt_clrp = "uc")


# plot results for t269
plotViolinplot(
  object = object,
  variables = genesetS2,
  across = "Lloyd_k7", 
  clrp = "npg",
  nrow = 3 
) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  legendTop()

############################## EXTRACT data frame of ssGSEA ###########################

UKF1_df <- Spata_GS[["data"]]
UKF2_df <- Spata_GS[["data"]]
UKF3_df <- Spata_GS[["data"]]
UKF4_df <- Spata_GS[["data"]]
UKF5_df <- Spata_GS[["data"]]
UKF6_df <- Spata_GS[["data"]]


############################### HEATMAP ######################################
library(dplyr)

# Add responder status to each filtered dataset
UKF1_df$patient <- "UKF1"
UKF2_df$patient <- "UKF2"
UKF3_df$patient <- "UKF3"
UKF4_df$patient <- "UKF4"
UKF5_df$patient <- "UKF5"
UKF6_df$patient <- "UKF6"

# Tag responder status
UKF1_df$responder_status <- "Non-Responder"
UKF2_df$responder_status <- "Non-Responder"
UKF3_df$responder_status <- "Non-Responder"
UKF4_df$responder_status <- "Responder"
UKF5_df$responder_status <- "Responder"
UKF6_df$responder_status <- "Responder"

# Combine into one data frame
combined_df <- bind_rows(
  UKF1_df, UKF2_df, UKF3_df,
  UKF4_df, UKF5_df, UKF6_df
)

###############################################################################
library(dplyr)
library(tidyr)
library(pheatmap)

# Step 1: average expression per variable × patient × responder_status
heatmap_df <- combined_df %>%
  group_by(responder_status, variables, patient) %>%
  summarise(mean_value = mean(values, na.rm = TRUE), .groups = "drop")

# Filter non-responders
df_nr <- heatmap_df %>%
  filter(responder_status == "Non-Responder") %>%
  select(variables, patient, mean_value) %>%
  pivot_wider(names_from = patient, values_from = mean_value) %>%
  column_to_rownames("variables") %>%
  as.matrix()

# Row-scale and impute low values
df_nr_scaled <- t(scale(t(df_nr)))
min_val_nr <- min(df_nr_scaled, na.rm = TRUE)
df_nr_scaled[is.na(df_nr_scaled)] <- min_val_nr / 100

# Plot
pheatmap(
  df_nr_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Non-Responders (UKF1–UKF3)",
  color = colorRampPalette(c("blue", "white", "red"))(100)
)



# Filter responders
df_r <- heatmap_df %>%
  filter(responder_status == "Responder") %>%
  select(variables, patient, mean_value) %>%
  pivot_wider(names_from = patient, values_from = mean_value) %>%
  column_to_rownames("variables") %>%
  as.matrix()

# Row-scale and impute low values
df_r_scaled <- t(scale(t(df_r)))
min_val_r <- min(df_r_scaled, na.rm = TRUE)
df_r_scaled[is.na(df_r_scaled)] <- min_val_r / 100

# Plot
pheatmap(
  df_r_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Responders (UKF4–UKF6)",
  color = colorRampPalette(c("blue", "white", "red"))(100)
)


########################### DeCouple R #########################################

library(readxl)
library(openxlsx)
library(dplyr)
library(tidyr)
library(decoupleR)

# Read Excel file (assumes sheet 1 and headers in first row)
lam_meta <- read_excel('/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATIAL_METABOLISM/LAM_Meta.xlsx')

# Convert wide format to long format (decoupleR net format)
net <- lam_meta %>%
  pivot_longer(cols = everything(), names_to = "source", values_to = "target") %>%
  filter(!is.na(target)) %>%
  mutate(weight = 1)

# View the result
head(net)

# Optional: save to CSV for decoupleR usage
write.xlsx(lam_net, "/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATIAL_METABOLISM/LAM_NET.xlsx")

# get matrix
# Check available matrix names
getMatrixNames(object)

# Get matrix — make sure to match an available name (e.g., "CLR", "counts")
mat <- getMatrix(object, mtr_name = "CLR", assay_name = activeAssay(object))

# Run mlm
acts <- decoupleR::run_ulm(mat = mat, 
                           net = net, 
                           .source = 'source', 
                           .target = 'target',
                           .mor = 'weight', 
                           minsize = 5)

#################################### add to meta data ################################################
meta <- object@meta_obs  # or: meta <- getCellData(object)

ulm_wide <- acts %>%
  select(condition, source, score) %>%
  pivot_wider(names_from = source, values_from = score)

# Ensure barcodes are rownames
ulm_meta <- ulm_wide %>%
  rename(barcodes = condition) %>%
  column_to_rownames("barcodes")

# Keep barcodes as a column, do NOT set as rownames
ulm_meta <- ulm_wide %>%
  rename(barcodes = condition)

object <- addMetaDataObs(
  object = object,
  meta_obs_df = ulm_meta,
  var_names = NULL,       # NULL = add all columns
  na_warn = TRUE,
  overwrite = FALSE       # Set to TRUE if replacing old scores
)

plotSurface(object, 
            color_by = "Lipid_Uptake_and_Transport", 
            smooth = FALSE, 
            smooth_span = 0.3)

# ssGSEA
Spata_GS <- plotSurfaceComparison(object = object, 
                                  color_by = "Lipid_Uptake_and_Transport",
                                  method_gs = "zscore",
                                  smooth = TRUE,
                                  pt_clrsp = "milo",
                                  smooth_span = 0.2,
                                  pt_size = 2,
                                  display_image = FALSE
)

Spata_GS

####################################################################################################
# Convert activity scores to wide format
acts_wide <- acts %>%
  select(condition, source, score) %>%
  pivot_wider(names_from = source, values_from = score)

# Merge ULM scores with metadata
merged_meta <- meta %>%
  left_join(ulm_wide, by = c("barcodes" = "condition"))

object <- setCellData(object, cell_data = merged_meta)


# acts_wide$condition == meta$barcodes
merged_df <- meta %>%
  left_join(acts_wide, by = c("barcodes" = "condition"))

write.xlsx(merged_df, '/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATIAL_METABOLISM/DECOUPLER/SPATIAL/UKF6.xlsx')

######################## HEATMAP ######################################
library(dplyr)
library(tidyr)
library(Seurat)
library(pheatmap)
library(RColorBrewer)
library(tibble)

selected_pathways <- c(
  "Cholesterol_Metabolism",
  "Fatty_Acid_Metabolism_General",
  "Fatty_Acid_Oxidation",
  "Fatty_Acid_Synthesis",
  "Glycolysis",
  "Lipid_Uptake_and_Transport",
  "Lipolysis",
  "Oxidative_Phosphorylation",
  "ROS_Production"
)


# Group by cluster and summarize average gene set scores
df <- merged_df %>%
  select(bayes_space, all_of(selected_pathways)) %>%
  pivot_longer(cols = -bayes_space, names_to = "source", values_to = "score") %>%
  group_by(bayes_space, source) %>%
  summarise(mean = mean(score, na.rm = TRUE), .groups = "drop")

# Pivot to wide matrix: clusters as rows, gene sets as columns
top_acts_mat <- df %>%
  pivot_wider(id_cols = bayes_space, names_from = source, values_from = mean) %>%
  column_to_rownames(var = "bayes_space") %>%
  as.matrix()

# Diverging RdBu palette
colors <- rev(brewer.pal(n = 11, name = "RdBu"))
colors.use <- colorRampPalette(colors)(100)

my_breaks <- c(seq(-1.25, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.05, 1.25, length.out = floor(100 / 2)))

pheatmap(mat = top_acts_mat,
         color = colors.use,
         breaks = my_breaks,
         border_color = "white",
         cellwidth = 20,
         cellheight = 20,
         treeheight_row = 20,
         treeheight_col = 20,
         cluster_rows = TRUE,
         cluster_cols = TRUE)

################################## TAM and Tumor HEATMAP ########################
library(readxl)
library(purrr)
library(stringr)

# Define path
folder_path <- "/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATIAL_METABOLISM/DECOUPLER/SPATIAL"

# List all Excel files (.xlsx or .xls)
file_list <- list.files(path = folder_path, pattern = "\\.xlsx?$", full.names = TRUE)

# Read all Excel files into a named list
excel_data_list <- file_list %>%
  set_names(~str_remove(basename(.x), "\\.xlsx?$")) %>%  # use file names as list names
  map(read_excel)

# Example: View one table
names(excel_data_list)  # List of datasets
head(excel_data_list[[1]])  # First Excel file's content
# Update sample column in each element of the list
excel_data_list <- imap(excel_data_list, ~ {
  .x$sample <- .y  # .y is the list name (e.g., "UKF1", "UKF2", ...)
  .x
})

library(dplyr)
library(purrr)

# Define cluster filters and labels
library(dplyr)
library(purrr)

# Define cluster filters and new region labels
cluster_filters <- list(
  UKF1 = list(groups = c(5, 1), labels = c(`5` = "LA-TAM", `1` = "Tumor")),
  UKF2 = list(groups = c(5, 4), labels = c(`5` = "LA-TAM", `4` = "Tumor")),
  UKF3 = list(groups = c(5, 1), labels = c(`5` = "LA-TAM", `1` = "Tumor")),
  UKF4 = list(groups = c(5, 7), labels = c(`5` = "IFN-TAM", `7` = "Tumor")),
  UKF5 = list(groups = c(10, 9), labels = c(`10` = "IFN-TAM", `9` = "Tumor")),
  UKF6 = list(groups = c(1, 6, 2), labels = c(`1` = "IFN-TAM", `6` = "IFN-TAM", `2` = "Tumor"))
)

# Apply filtering and assign 'region' based on cluster
filtered_list <- imap(excel_data_list, function(df, sample_name) {
  filter_info <- cluster_filters[[sample_name]]
  df %>%
    filter(Cluster %in% filter_info$groups) %>%
    mutate(region = recode(as.character(Cluster), !!!filter_info$labels))
})

# heatmap
library(dplyr)
library(tidyr)
library(purrr)

# Combine filtered list into a single data frame
# Coerce Cluster column to character in each element
filtered_list <- map(filtered_list, ~ .x %>%
                       mutate(Cluster = as.character(Cluster)))
merged_filtered_df <- bind_rows(filtered_list, .id = "sample_name")

# Create a unique row label: Sample-Region
merged_filtered_df <- merged_filtered_df %>%
  mutate(sample_region = paste(sample_name, region, sep = "-"))

# List of pathways (as per your earlier column names)
pathways <- c(
  "Cholesterol_Metabolism",
  "Fatty_Acid_Oxidation",
  "Fatty_Acid_Synthesis",
  "Glycolysis",
  "Lipid_Uptake_and_Transport",
  "Lipolysis",
  "Oxidative_Phosphorylation",
  "ROS_Production"
)

# Compute average scores per sample-region
heatmap_df <- merged_filtered_df %>%
  group_by(sample_region) %>%
  summarise(across(all_of(pathways), ~mean(.x, na.rm = TRUE))) %>%
  ungroup() %>%
  column_to_rownames("sample_region")


write.xlsx(heatmap_df, "/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATIAL_METABOLISM/DECOUPLER/SPATIAL/heatmap_df.xlsx", rowNames = TRUE)

#if need be
heatmap_df <- read.xlsx("/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATIAL_METABOLISM/DECOUPLER/SPATIAL/heatmap_df.xlsx", rowNames = TRUE)

# Scale each column (z-score), then replace NAs with 0
heatmap_scaled <- scale(heatmap_df)
heatmap_scaled[is.na(heatmap_scaled)] <- 0

library(pheatmap)
library(RColorBrewer)

# Color palette
colors <- rev(brewer.pal(n = 11, name = "RdBu"))
colors.use <- colorRampPalette(colors)(100)

# Plot heatmap
pheatmap(mat = heatmap_scaled,
         color = colors.use,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         border_color = "white",
         cellwidth = 25,
         cellheight = 20,
         treeheight_row = 0,
         treeheight_col = 20,
         fontsize = 10)

############################### GOI surface plot ####################################

Lipid_up <- c("LDLR", "LRP1", "CD36", "MSR1", "SCARB1")

Lipid_up_UKF3 <- c("CD36", "FABP4", "FABP5")

Lipolysis <- c("PNPLA2", "LIPE", "MGLL", "LPL", "LIPA")

CholMeta <- c("HMGCR",	"SQLE",	"DHCR24",	"ACAT1",	"ACAT2",	"ABCA1",	"ABCG1",	"NR1H2",	"NR1H3",	"SREBF1",	"SREBF2",	"NPC1",	"NPC2")

# ssGSEA
Spata_GS <- plotSurfaceComparison(object = object, 
                                  color_by = Lipid_up_UKF3,
                                  method_gs = "zscore",
                                  smooth = TRUE,
                                  smooth_span = 0.2,
                                  pt_size = 3,
                                  display_image = FALSE
)

Spata_GS

############################################# PATHWAY MAPPPING #######################################################
BiocManager::install(c("pathview", "KEGGREST", "org.Hs.eg.db"))
if (!requireNamespace("gage", quietly = TRUE))
  install.packages("gage")
if (!requireNamespace("gageData", quietly = TRUE))
  BiocManager::install("gageData")

library(pathview)
library(KEGGREST)

genes <- c("LDLR",
           "LRP1",
           "CD36",
           "HMGCR",
           "SQLE",
           "DHCR24",
           "ACAT1",
           "PNPLA2",
           "MGLL",
           "LIPA",
           "LPL",
           "APOE", 
           "TREM2")

# Get matrix — make sure to match an available name (e.g., "CLR", "counts")
mat <- getMatrix(object, mtr_name = "CLR", assay_name = activeAssay(object))

# Ensure gene names are matched exactly with rownames
genes_to_use <- intersect(genes, rownames(mat))

# Subset the matrix
mat_subset <- mat[genes_to_use, , drop = FALSE]

meta <- object@meta_obs

# Create groupings based on bayes_space
group_c5 <- meta$barcodes[meta$bayes_space == 5]
group_rest <- meta$barcodes[meta$bayes_space != 5]

# Make sure barcodes match columns of the matrix
group_c5 <- intersect(group_c5, colnames(mat_subset))
group_rest <- intersect(group_rest, colnames(mat_subset))

# Calculate average expression per gene for each group
avg_c5 <- rowMeans(mat_subset[, group_c5, drop = FALSE])
avg_rest <- rowMeans(mat_subset[, group_rest, drop = FALSE])

# Compute log2 fold change (add pseudocount to avoid log(0))
log2_fc <- log2((avg_c5 + 1e-6) / (avg_rest + 1e-6))

# Create a data frame with results
fc_results <- data.frame(
  gene = rownames(mat_subset),
  avg_C5 = avg_c5,
  avg_rest = avg_rest,
  log2FC = log2_fc
)

# Optional: sort by log2FC
fc_results <- fc_results[order(-fc_results$log2FC), ]

################################## MAP #######################################

# Assuming fc_results is a data.frame with rownames as gene symbols
log2fc_vector <- fc_results$log2FC
names(log2fc_vector) <- rownames(fc_results)

library(org.Hs.eg.db)

entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = names(log2fc_vector),
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

# Clean up: remove NAs
log2fc_entrez <- log2fc_vector[!is.na(entrez_ids)]
names(log2fc_entrez) <- entrez_ids[!is.na(entrez_ids)]

library(KEGGREST)
kegg_pathways <- keggList("pathway", "hsa")
grep("lipid|cholesterol|fatty", kegg_pathways, ignore.case = TRUE, value = TRUE)



# Example: convert gene symbols to Entrez IDs
library(org.Hs.eg.db)
genes <- rownames(fc_results)
entrez_ids <- mapIds(org.Hs.eg.db, keys = genes,
                     column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

# Remove NA mappings
entrez_ids <- entrez_ids[!is.na(entrez_ids)]

# Query KEGG for each gene's pathways
gene2pathway <- lapply(entrez_ids, function(eid) {
  keggLink("pathway", paste0("hsa:", eid))
})

# Organize as a data.frame
pathway_map <- data.frame(
  Gene = rep(names(gene2pathway), lengths(gene2pathway)),
  Entrez = rep(entrez_ids, lengths(gene2pathway)),
  Pathway = unlist(gene2pathway)
)

head(pathway_map)

############################################ VISUALIAZION ###################
#BiocManager::install(c("rWikiPathways", "RCy3"))
library(rWikiPathways)
library(RCy3)

# Search for cholesterol metabolism pathway
results <- findPathwaysByText("cholesterol metabolism")
# Filter for Homo sapiens
results <- results[results$species == "Homo sapiens", ]
results

# Make sure Cytoscape is running before this script
cytoscapePing()

setwd('/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATIAL_METABOLISM/pathway_visualization')

# Download WP5304 in GPML format
gpml <- getPathway("WP5304")
writeLines(gpml, "WP5304.gpml")

library(xml2)
library(dplyr)
library(ggplot2)

# Parse GPML
doc <- read_xml(gpml)

# Extract DataNodes
nodes <- xml_find_all(doc, ".//d1:DataNode", xml_ns(doc))

# Extract info from both the DataNode and its Graphics child
node_df <- data.frame(
  label = xml_attr(nodes, "TextLabel"),
  x = as.numeric(xml_attr(xml_find_first(nodes, ".//d1:Graphics", xml_ns(doc)), "CenterX")),
  y = as.numeric(xml_attr(xml_find_first(nodes, ".//d1:Graphics", xml_ns(doc)), "CenterY")),
  stringsAsFactors = FALSE
)

# Join FC values
node_df <- node_df %>%
  left_join(
    fc_results %>% dplyr::select(gene, log2FC),
    by = c("label" = "gene")
  )

# Plot with coloring
ggplot(node_df, aes(x = x, y = -y, label = label, fill = log2FC)) +
  geom_point(shape = 21, size = 5) +
  geom_text(vjust = -1, size = 3) +
  scale_fill_gradient(low = "#FFCCCC", high = "#990000", na.value = "grey90") +
  theme_void() +
  labs(title = "WP5304 - Cholesterol Pathway (log2FC coloring)",
       fill = "log2FC")

######################## Enrichment score violin plot NonRes vs. Res #####################################

# load data 
data <- read.xlsx('/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATIAL_METABOLISM/NonRes_vs_Res_mata_Enrichment_values.xlsx')

# load libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# Reshape to long format
data_long <- data %>%
  pivot_longer(cols = c(LATAM, Tumor),
               names_to = "Compartment",
               values_to = "Score")

# Plot: side-by-side violins for LATAM vs Tumor, grouped by X1
ggplot(data_long, aes(x = Compartment, y = Score, fill = Status)) +
  geom_violin(trim = FALSE, position = position_dodge(width = 0.9), alpha = 0.7, width = 0.8) +
  geom_jitter(aes(color = Status),
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9),
              size = 1.8, alpha = 0.8) +
  geom_boxplot(width = 0.1, outlier.shape = NA,
               position = position_dodge(width = 0.9), alpha = 0.4) +
  stat_compare_means(aes(group = Status),
                     method = "t.test",
                     label = "p.signif",
                     label.y.npc = 0.95,
                     size = 3) +
  facet_wrap(~ X1, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Responder" = "#4CAF50", "NonResponder" = "#F44336")) +
  scale_color_manual(values = c("Responder" = "#388E3C", "NonResponder" = "#D32F2F")) +
  theme_bw(base_size = 12) +
  labs(title = "Side-by-side comparison of LATAM vs Tumor per metabolic term",
       x = "Compartment", y = "Score") +
  theme(
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "gray95"),
    panel.grid.major = element_line(color = "gray90"),
    legend.position = "top"
  )

################################### DEcpupleR ULM enrichment REs vs. NonREs#######################
# load libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(openxlsx)

data <- read.xlsx('/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATIAL_METABOLISM/DECOUPLER/SPATIAL/heatmap_df.xlsx')

library(tidyverse)
library(ggpubr)

data_long <- data %>%
  separate(X1, into = c("Sample", "CellType"), sep = "-", remove = FALSE) %>%
  mutate(
    Response = case_when(
      Sample %in% c("UKF1", "UKF2", "UKF3") ~ "NonResponder",
      Sample %in% c("UKF4", "UKF5", "UKF6") ~ "Responder",
      TRUE ~ NA_character_
    )
  ) %>%
  pivot_longer(
    cols = -c(X1, Sample, CellType, Response),
    names_to = "Pathway",
    values_to = "Score"
  )

ggplot(data_long, aes(x = Response, y = Score, fill = Response)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 1) +
  stat_summary(fun = "median", geom = "point", size = 2, color = "black") +
  stat_compare_means(
    aes(group = Response),
    method = "wilcox.test",
    label = "p.signif"
  ) +
  facet_grid(CellType ~ Pathway, scales = "free_y") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text.x = element_text(size = 8),
    strip.text.y = element_text(size = 9)
  )


################################# DecoupleR FIANL #####################################
library(dplyr)
library(tidyr)
library(stringr)

# load net file 
net <- read.xlsx('/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATIAL_METABOLISM/LAM_NET.xlsx')
# load cluster info
cluster <- read.xlsx('/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATIAL_METABOLISM/Cluster_Annotation.xlsx')

# Ensure cluster_df (the annotation file) is loaded
# Columns: Sample | Annotation | ClusterA | ClusterB
cluster_df <- cluster_df %>%
  tidyr::pivot_longer(cols = starts_with("Cluster"),
                      names_to = "ClusterType",
                      values_to = "Cluster") %>%
  filter(!is.na(Cluster)) %>%
  select(Sample, Annotation, Cluster) %>%
  mutate(Cluster = as.character(Cluster))

# Define cluster column names per sample
cluster_columns <- list(
  UKF1 = "bayes_space",
  UKF2 = "bayes_space",
  UKF3 = "Lloyd_k6",
  UKF4 = "bayes_space",
  UKF5 = "bayes_space",
  UKF6 = "Lloyd_k7"
)

# --- Function to rename cluster IDs to biological annotation ---
rename_clusters_in_spata <- function(spata_obj, sample_name, cluster_df, cluster_columns){
  
  # get mapping of cluster IDs → Annotation for this sample
  mapping_df <- cluster_df %>%
    filter(Sample == sample_name) %>%
    mutate(Cluster = as.character(Cluster))
  
  # detect which metadata column holds cluster assignments
  cluster_col <- cluster_columns[[sample_name]]
  meta_df <- getMetaDf(spata_obj)
  
  if(!cluster_col %in% colnames(meta_df)){
    stop(paste0("❌ Column ", cluster_col, " not found in ", sample_name))
  }
  
  # clean cluster labels (strip letters)
  meta_df <- meta_df %>%
    mutate(cluster_clean = str_remove(!!sym(cluster_col), "^[A-Za-z]"))
  
  # join to annotation
  meta_df <- meta_df %>%
    left_join(mapping_df, by = c("cluster_clean" = "Cluster")) %>%
    mutate(cluster_label = ifelse(!is.na(Annotation), Annotation, cluster_clean))
  
  # Update SPATA object metadata
  spata_obj@meta_obs$cluster_label <- meta_df$cluster_label
  
  message(paste0("✅ Renamed clusters in ", sample_name,
                 " (using column '", cluster_col, "')"))
  
  return(spata_obj)
}

# --- Apply to all SPATA subset objects ---
spata_subset_renamed <- purrr::imap(spata_subset_list,
                                    ~rename_clusters_in_spata(.x, .y, cluster_df, cluster_columns))


# Ensure your list is correctly ordered and named
names(spata_subset_renamed) <- c("UKF1", "UKF2", "UKF3", "UKF4", "UKF5", "UKF6")

# Rename the 'sample' metadata field and internal slot
spata_subset_renamed <- imap(spata_subset_renamed, function(spata_obj, sample_name) {
  
  # Update sample column in metadata
  if("sample" %in% colnames(spata_obj@meta_obs)){
    spata_obj@meta_obs$sample <- sample_name
  } else {
    spata_obj@meta_obs <- spata_obj@meta_obs %>%
      mutate(sample = sample_name)
  }
  
  # Update internal sample slot (SPATA2 uses this for labeling)
  spata_obj@sample <- sample_name
  
  message(paste0("✅ Updated sample name to ", sample_name))
  return(spata_obj)
})

# Extract matrix
# Get matrix — make sure to match an available name (e.g., "CLR", "counts")
UKF1 <- spata_subset_renamed$UKF1
UKF2 <- spata_subset_renamed$UKF2
UKF3 <- spata_subset_renamed$UKF3
UKF4 <- spata_subset_renamed$UKF4
UKF5 <- spata_subset_renamed$UKF5
UKF6 <- spata_subset_renamed$UKF6


mat1 <- getMatrix(UKF1, mtr_name = "CLR", assay_name = activeAssay(UKF1))
mat2 <- getMatrix(UKF2, mtr_name = "CLR", assay_name = activeAssay(UKF2))
mat3 <- getMatrix(UKF3, mtr_name = "CLR", assay_name = activeAssay(UKF3))
mat4 <- getMatrix(UKF4, mtr_name = "CLR", assay_name = activeAssay(UKF4))
mat5 <- getMatrix(UKF5, mtr_name = "CLR", assay_name = activeAssay(UKF5))
mat6 <- getMatrix(UKF6, mtr_name = "CLR", assay_name = activeAssay(UKF6))

# Run mlm
acts6 <- decoupleR::run_ulm(mat = mat6, 
                           net = net, 
                           .source = 'source', 
                           .target = 'target',
                           .mor = 'weight', 
                           minsize = 1)

# addmeta data to acts df
library(dplyr)
library(purrr)

# Combine all your acts tables
acts_list <- list(
  UKF1 = acts1,
  UKF2 = acts2,
  UKF3 = acts3,
  UKF4 = acts4,
  UKF5 = acts5,
  UKF6 = acts6
)

# Function to attach sample and cluster_label
add_sample_and_cluster <- function(acts_tbl, spata_obj, sample_name) {
  
  # Ensure consistent column naming
  acts_tbl <- acts_tbl %>%
    rename(barcodes = condition)
  
  # Extract barcode–annotation mapping
  meta_info <- spata_obj@meta_obs %>%
    dplyr::select(barcodes, cluster_label)
  
  # Merge metadata
  acts_tbl <- acts_tbl %>%
    left_join(meta_info, by = "barcodes") %>%
    mutate(sample = sample_name)
  
  message(paste0("✅ Added sample + cluster_label to ", sample_name,
                 " (", nrow(acts_tbl), " entries)"))
  
  return(acts_tbl)
}

# Apply to all samples
acts_annotated_list <- purrr::imap(acts_list, function(acts_tbl, sample_name){
  add_sample_and_cluster(acts_tbl, spata_subset_renamed[[sample_name]], sample_name)
})

# Optionally combine into one big data frame
acts_combined <- dplyr::bind_rows(acts_annotated_list)

# Inspect
head(acts_combined)

write.xlsx(acts_combined, "/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATIAL_METABOLISM/DECOUPLER/SPATIAL/acts_combined_df.xlsx")


# Assuming acts_combined looks like:
# statistic | source | barcodes | score | p_value | cluster_label | sample

# 1️⃣ Split the full table into a list of 6 data frames by sample
acts_split <- split(acts_combined, acts_combined$sample)

# 2️⃣ For each sample, pivot wider so 'source' becomes column names
acts_wide_list <- purrr::imap(acts_split, function(df, sample_name) {
  
  df_wide <- df %>%
    select(barcodes, source, score, cluster_label, sample) %>%
    pivot_wider(
      names_from  = source,   # each pathway becomes a column
      values_from = score
    ) %>%
    distinct(barcodes, .keep_all = TRUE)
  
  message(paste0("✅ Created wide-format ULM data for ", sample_name, 
                 " with ", nrow(df_wide), " barcodes and ", ncol(df_wide)-3, " pathways."))
  
  return(df_wide)
})

# 3️⃣ Access individual data frames
acts_UKF1 <- acts_wide_list$UKF1
acts_UKF2 <- acts_wide_list$UKF2
acts_UKF3 <- acts_wide_list$UKF3
acts_UKF4 <- acts_wide_list$UKF4
acts_UKF5 <- acts_wide_list$UKF5
acts_UKF6 <- acts_wide_list$UKF6

# (Optional) check one
head(acts_UKF1[, 1:8])


# mod dfs




library(openxlsx)
library(purrr)

# Define output folder
out_dir <- "/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATIAL_METABOLISM/DECOUPLER/SPATIAL/FINAL_scts_dfs"

# Make sure folder exists
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Combine all six data frames into a named list
acts_wide_list <- list(
  UKF1 = acts_UKF1,
  UKF2 = acts_UKF2,
  UKF3 = acts_UKF3,
  UKF4 = acts_UKF4,
  UKF5 = acts_UKF5,
  UKF6 = acts_UKF6
)

# Write each one to its own Excel file (simple format)
purrr::iwalk(acts_wide_list, function(df, name) {
  file_path <- file.path(out_dir, paste0(name, "_ULM_scores.xlsx"))
  writexl::write_xlsx(df, path = file_path)
  message(paste0("✅ Wrote ", file_path))
})




################################ Genes comparsion LA-TAM vs. Tumor ####################
# ======================================================
# 1. Load libraries
# ======================================================
library(SPATA2)
library(dplyr)
library(purrr)
library(tidyr)
library(openxlsx)

# ======================================================
# 2. Define gene list of interest
# ======================================================
genes <- c("TREM2", "LDLR", "LRP1", "CD36", "PNPLA2",
           "MGLL", "LPL", "ACAT1", "LIPA", "DHCR24",
           "APOE", "HMGCR", "SQLE")

# ======================================================
# 3. Load SPATA objects (replace with your actual objects)
# ======================================================
spata_list <- list(
  UKF1 = UKF1,
  UKF2 = UKF2,
  UKF3 = UKF3,
  UKF4 = UKF4,
  UKF5 = UKF5,
  UKF6 = UKF6
)

# ======================================================
# 4. Load and reshape cluster annotation
# ======================================================
cluster_df <- read.xlsx("/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATIAL_METABOLISM/Cluster_Annotation.xlsx")

# Combine ClusterA and ClusterB into one column
cluster_long <- cluster_df %>%
  pivot_longer(cols = starts_with("Cluster"),
               names_to = "ClusterType",
               values_to = "Cluster") %>%
  filter(!is.na(Cluster)) %>%
  select(Sample, Annotation, Cluster)

# ======================================================
# 5. Define function to extract gene expression + metadata
# ======================================================
extract_gene_data <- function(obj_name, spata_obj, cluster_long, genes){
  
  # 1. Detect available matrices and pick the first valid one
  mtr_names <- getMatrixNames(spata_obj)
  preferred_order <- c("CLR", "scaled", "normalized", "logcounts", "counts")
  mtr_name <- intersect(preferred_order, mtr_names)[1]
  
  if(is.na(mtr_name)){
    stop(paste0("❌ No valid matrix found in ", obj_name, ". Available: ", paste(mtr_names, collapse = ", ")))
  } else {
    message(paste0("✅ Using matrix '", mtr_name, "' for ", obj_name))
  }
  
  # 2. Detect which genes are available
  available_genes <- getGenes(spata_obj)
  genes_present <- intersect(genes, available_genes)
  
  if(length(genes_present) == 0){
    message(paste0("⚠️ No target genes found in ", obj_name, ". Skipping."))
    return(NULL)
  }
  
  # 3. Extract the selected expression matrix
  expr_mtr <- getMatrix(spata_obj, mtr_name = mtr_name)
  expr_df <- expr_mtr[genes_present, , drop = FALSE] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(-gene, names_to = "barcodes", values_to = "expression")
  
  # 4. Extract metadata
  meta_df <- getMetaDf(spata_obj)
  
  # 5. Detect correct cluster column (varies across samples)
  cluster_col <- case_when(
    "bayes_space" %in% colnames(meta_df) ~ "bayes_space",
    "Lloyd_k6" %in% colnames(meta_df) ~ "Lloyd_k6",
    "Lloyd_k7" %in% colnames(meta_df) ~ "Lloyd_k7",
    TRUE ~ NA_character_
  )
  
  if (is.na(cluster_col)) {
    stop(paste0("❌ No recognized cluster column found for ", obj_name))
  }
  
  meta_df <- meta_df %>%
    dplyr::select(barcodes, cluster = all_of(cluster_col)) %>%
    dplyr::mutate(cluster = as.character(cluster))
  
  # 6. Merge expression + metadata
  df <- expr_df %>%
    dplyr::left_join(meta_df, by = "barcodes") %>%
    dplyr::mutate(Sample = obj_name)
  
  # 7. Convert cluster_long$Cluster to character for type-safe joining
  cluster_long <- cluster_long %>%
    dplyr::mutate(Cluster = as.character(Cluster))
  
  # 8. Join with cluster annotation (LA-TAM, IFN-TAM, Tumor)
  df <- df %>%
    dplyr::left_join(cluster_long, by = c("Sample" = "Sample", "cluster" = "Cluster"))
  
  return(df)
}



# ======================================================
# 6. Apply to all SPATA objects
# ======================================================
combined_df <- purrr::imap_dfr(spata_list, ~extract_gene_data(.y, .x, cluster_long, genes))

# ======================================================
# 7. Filter relevant annotations and summarize
# ======================================================
subset_df <- combined_df %>%
  filter(Annotation %in% c("LA-TAM", "IFN-TAM", "Tumor"))

summary_df <- combined_df %>%
  filter(Annotation %in% c("LA-TAM", "IFN-TAM", "Tumor")) %>%
  group_by(Sample, Annotation, gene) %>%
  summarise(mean_expr = mean(expression, na.rm = TRUE), .groups = "drop")


# ======================================================
# 8. Inspect results
# ======================================================
head(summary_df)


summary_df <- combined_df %>%
  filter(Annotation %in% c("LA-TAM", "IFN-TAM", "Tumor")) %>%
  group_by(Sample, Annotation, gene) %>%
  summarise(mean_expr = mean(expression, na.rm = TRUE), .groups = "drop")

# Pivot wider so annotations become columns
fc_df <- summary_df %>%
  tidyr::pivot_wider(
    names_from = Annotation,
    values_from = mean_expr
  ) %>%
  # Compute fold changes per sample
  mutate(
    FC_LA_TAM_vs_Tumor = `LA-TAM` / Tumor,
    FC_IFN_TAM_vs_Tumor = `IFN-TAM` / Tumor
  )

# Average FC across all samples for each gene
avg_fc_df <- fc_df %>%
  group_by(gene) %>%
  summarise(
    avg_FC_LA_TAM_vs_Tumor = mean(FC_LA_TAM_vs_Tumor, na.rm = TRUE),
    avg_FC_IFN_TAM_vs_Tumor = mean(FC_IFN_TAM_vs_Tumor, na.rm = TRUE)
  ) %>%
  ungroup()

# Inspect result
avg_fc_df
write.xlsx(avg_fc_df, "/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATIAL_METABOLISM/LATAM_IFNTAM_Tumor_FC.xlsx")


################################ FINAL LATAM vs. IFNTAM vs. Tumor ###########################################
# laod excel files of ULM values which are only significant
library(readxl)
library(dplyr)
library(purrr)

# Define your folder path
path <- "/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/SpatIal_analysis_02.02.2024/SPATIAL_METABOLISM/DECOUPLER/SPATIAL/FINAL_scts_dfs"

# List all Excel files (.xlsx)
files <- list.files(path, pattern = "\\.xlsx$", full.names = TRUE)

# Read all files into a named list
acts_list <- purrr::map(files, readxl::read_xlsx) %>%
  purrr::set_names(basename(files))  # keeps filenames as list names

# Check
names(acts_list)
length(acts_list)

acts_combined <- purrr::imap_dfr(acts_list, ~mutate(.x, sample = tools::file_path_sans_ext(.y)))

# Inspect
dim(acts_combined)
head(acts_combined[, 1:8])


# plot and stat
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# --- Clean sample names ---
acts_combined <- acts_combined %>%
  mutate(sample = gsub("_ULM_scores", "", sample))

# --- Define Non-Responders and Responders ---
nonres_samples <- c("UKF1", "UKF2", "UKF3")
res_samples    <- c("UKF4", "UKF5", "UKF6")

# --- Pathway columns ---
pathway_cols <- c(
  "Cholesterol_Metabolism",
  "Fatty_Acid_Oxidation",
  "Fatty_Acid_Synthesis",
  "Glycolysis",
  "Lipid_Uptake_and_Transport",
  "Lipolysis",
  "Oxidative_Phosphorylation",
  "ROS_Production"
)

# --- Reshape to long format ---
df_long <- acts_combined %>%
  select(sample, cluster_label, all_of(pathway_cols)) %>%
  pivot_longer(cols = all_of(pathway_cols),
               names_to = "Pathway",
               values_to = "Score")

# --- Assign biological group labels ---
df_long <- df_long %>%
  mutate(
    Cohort = case_when(
      sample %in% nonres_samples ~ "Non-Responder",
      sample %in% res_samples ~ "Responder",
      TRUE ~ NA_character_
    ),
    Group = case_when(
      Cohort == "Non-Responder" & cluster_label == "LA-TAM" ~ "LA-TAM",
      Cohort == "Non-Responder" & cluster_label == "Tumor"  ~ "Tumor_NonRes",
      Cohort == "Responder" & cluster_label == "IFN-TAM"   ~ "IFN-TAM",
      Cohort == "Responder" & cluster_label == "Tumor"     ~ "Tumor_Res",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group))

# --- Set order for x-axis ---
df_long$Group <- factor(df_long$Group,
                        levels = c("LA-TAM", "Tumor_NonRes", "IFN-TAM", "Tumor_Res"))

# --- Clean pathway labels for nice facets ---
df_long$Pathway <- gsub("_", " ", df_long$Pathway)

# --- Colors for each group ---
cols <- c(
  "LA-TAM" = "#FF4041",
  "Tumor_NonRes" = "#FAAB33",
  "IFN-TAM" = "#487BE3",
  "Tumor_Res" = "#F99244"
)

# --- Define all biologically meaningful comparisons ---
my_comparisons <- list(
  c("LA-TAM", "Tumor_NonRes"),
  c("IFN-TAM", "Tumor_Res"),
  c("LA-TAM", "IFN-TAM"),
  c("Tumor_NonRes", "Tumor_Res")
)

# --- Plot ---
ggplot(df_long, aes(x = Group, y = Score, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.1) +
  scale_fill_manual(values = cols) +
  facet_wrap(~Pathway, scales = "free_y", ncol = 2) +
  theme_classic(base_size = 13) +
  labs(x = NULL, y = "Pathway enrichment (spot-level)") +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.format",
    hide.ns = TRUE
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 10, angle = 15, hjust = 1),
    strip.text = element_text(face = "bold", size = 11)
  )

# ----- vionlin box plot ---------
library(ggplot2)
library(ggpubr)

# --- Define the desired facet order ---
df_long$Pathway <- factor(df_long$Pathway, levels = c(
  "Lipid Uptake and Transport",
  "Lipolysis",
  "Cholesterol Metabolism",
  "Fatty Acid Oxidation",
  "Fatty Acid Synthesis",
  "Oxidative Phosphorylation",
  "Glycolysis",
  "ROS Production"
))

# --- Plot ---
ggplot(df_long, aes(x = Group, y = Score, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.8, width = 1, color = "gray60") +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7, color = "black") +
  scale_fill_manual(values = cols) +
  facet_wrap(~Pathway, scales = "free_y", ncol = 4) +
  theme_classic(base_size = 13) +
  labs(x = NULL, y = "Pathway enrichment (spot-level)") +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.format",
    hide.ns = TRUE
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 10, angle = 15, hjust = 1),
    strip.text = element_text(face = "bold", size = 11),
    panel.spacing = unit(0.8, "lines")
  )


########################################### Pathway score per patient ##################################
library(dplyr)
library(tidyr)
library(ggplot2)

# --- Prepare data ---
df_long <- acts_combined %>%
  filter(cluster_label %in% c("LA-TAM", "IFN-TAM")) %>%
  select(sample, cluster_label,
         Cholesterol_Metabolism, Fatty_Acid_Oxidation, Fatty_Acid_Synthesis,
         Glycolysis, Lipid_Uptake_and_Transport, Lipolysis,
         Oxidative_Phosphorylation, ROS_Production) %>%
  pivot_longer(
    cols = Cholesterol_Metabolism:ROS_Production,
    names_to = "Pathway",
    values_to = "Score"
  ) %>%
  group_by(sample, cluster_label, Pathway) %>%
  summarise(
    Mean_Score = mean(Score, na.rm = TRUE),
    SEM = sd(Score, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# --- Set factor levels for consistent order ---
path_order <- c("Lipid_Uptake_and_Transport", "Lipolysis",
                "Cholesterol_Metabolism", "Fatty_Acid_Oxidation",
                "Fatty_Acid_Synthesis", "Oxidative_Phosphorylation",
                "Glycolysis", "ROS_Production")
df_long$Pathway <- factor(df_long$Pathway, levels = path_order)

# UKF1–3 (NR) at top, UKF4–6 (R) below
patient_order <- c("UKF1", "UKF2", "UKF3", "UKF4", "UKF5", "UKF6")
df_long$sample <- factor(df_long$sample, levels = patient_order)

# --- Plot colors ---
cols <- c("LA-TAM" = "#E74C3C", "IFN-TAM" = "#3498DB")

# --- Plot: unidirectional grouped barplot ---
ggplot(df_long, aes(x = sample, y = Mean_Score, fill = cluster_label)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(
    ymin = Mean_Score - SEM,
    ymax = Mean_Score + SEM
  ), 
  position = position_dodge(width = 0.7), width = 0.25, color = "black", size = 0.3) +
  facet_wrap(~ Pathway, scales = "free_x", ncol = 4) +
  scale_fill_manual(values = cols) +
  theme_classic(base_size = 13) +
  labs(
    x = "Mean pathway enrichment (± SEM)",
    y = NULL,
    fill = NULL
  ) +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 9),
    legend.position = "top"
  )


########################################### Pathway visualization ######################################
library(pathview)

# Vector of gene values (example: log2FC)
gene_data <- c(
  TREM2 = 1.2, LDLR = 2.1, LRP1 = 1.7, CD36 = 2.3, PNPLA2 = 0.8,
  MGLL = 1.5, LPL = 2.8, ACAT1 = 0.9, LIPA = 1.2, DHCR24 = 1.0,
  APOE = 2.4, HMGCR = 0.7, SQLE = 1.1
)

# Pathway ID for Cholesterol metabolism
pathview(
  gene.data = gene_data,
  pathway.id = "hsa04979",          # Cholesterol metabolism
  species = "hsa",                  # human
  limit = list(gene = 3),           # color scale limits
  low = list(gene = "white"), 
  mid = list(gene = "pink"), 
  high = list(gene = "red"),
  out.suffix = "custom_cholesterol"
)

