# load library
##############
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(patchwork)
library(vroom)
library(ggplot2)
library(openxlsx)

#####
############################################ ############################################ ####################################
##################################### Seurat cell type annotation in Tirosch data ############################################
############################################ ############################################ ####################################
#load
seurat_object <- readRDS('/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/Melnaoma_PBMC_scRNA/Tirosh_CAF_LA_TAM_cellchat_NichNet/seurat_object_celltypeFinal.rds')

#update
seurat_object <- UpdateSeuratObject(seurat_object)

# combine CAF
seurat_object@meta.data <- seurat_object@meta.data %>%
  mutate(
    celltypeFinal2 = case_when(
      celltypeFinal %in% c("contractile_CAF", "immune_CAF") ~ "CAF",
      celltypeFinal == "nonLA_TAM" ~ "IFN_TAM",
      TRUE ~ celltypeFinal
    )
  )

#set Idents
Idents(seurat_object) <- seurat_object$celltypeFinal2

custom_clolors <- c("Transitory" = "#ece30a",
                    "Neural.crest-like" = "#004f5e",
                    "Melanocytic" = "#009347",
                    "LA_TAM" = "#fca407",
                    "IFN_TAM" = "#4284e4",
                    "Exhausted_CD8" = "#f174da",
                    "Activated_CD8" = "#574073")


#umap
DimPlot(
  seurat_object,
  reduction = "umap",
  cols = custom_clolors,
  label = TRUE
)

#save
saveRDS(seurat_object, "/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/Melanoma_LATAM_Manuscript/Genome_Medicine_2026/Manuscript_Chakraborty_et_al_2026/Melanoma_Tirosch_cell_Annotation_CIBERSORTx/Tirosch_seurat_object_celltype.rds")

# extract inputs from seurat object for CIBERSORTx
Idents(seurat_object) <- seurat_object$celltypeFinal2
DefaultAssay(seurat_object) <- "RNA"
expr <- GetAssayData(
  seurat_object,
  assay = "RNA",
  layer = "data"
)
expr <- as.matrix(expr)
dim(expr)
head(rownames(expr))
head(colnames(expr))
cell_labels <- seurat_object$celltypeFinal2
names(cell_labels) <- colnames(seurat_object)
cell_labels <- cell_labels[colnames(expr)]
expr <- as.matrix(expr)
cibersort_mat <- cbind(GeneSymbol = rownames(expr), expr)
colnames(cibersort_mat)[-1] <- cell_labels

setwd('/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/Melanoma_LATAM_Manuscript/Genome_Medicine_2026/Manuscript_Chakraborty_et_al_2026/Melanoma_Tirosch_cell_Annotation_CIBERSORTx')

write.table(
  cibersort_mat,
  file = "CIBERSORTx_scRNA_reference.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


######################################### prepare mixture mat input for CIBERSORTX with meta data ############################
# sig_mat 
mixture_mat <- read.xlsx('/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/Melanoma_LATAM_Manuscript/Genome_Medicine_2026/Manuscript_Chakraborty_et_al_2026/Melanoma_Tirosch_cell_Annotation_CIBERSORTx/CIBERSORTX_input/Merged_Mat.xlsx')
mixture_meta <- read.xlsx('/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/Melanoma_LATAM_Manuscript/Genome_Medicine_2026/Manuscript_Chakraborty_et_al_2026/Melanoma_Tirosch_cell_Annotation_CIBERSORTx/CIBERSORTX_input/sig_mat_meta.xlsx')
# if mixture_mat is a data.frame with a Gene column
rownames(mixture_mat) <- mixture_mat$Gene

# remove the Gene column
mixture_mat$Gene <- NULL

# rename column
library(dplyr)

mixture_meta <- mixture_meta %>%
  rename(sample_id = Gene)

#subset m,atrix
mixture_mat <- mixture_mat[, colnames(mixture_mat) %in% mixture_meta$sample_id]
Mixture <- cbind(Gene = rownames(mixture_mat), mixture_mat)

# optional: remove rownames
rownames(Mixture) <- NULL

#save
write.xlsx(Mixture, "/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/Melanoma_LATAM_Manuscript/Genome_Medicine_2026/Manuscript_Chakraborty_et_al_2026/Melanoma_Tirosch_cell_Annotation_CIBERSORTx/CIBERSORTX_input/Mixture_Mat.xlsx")

setwd('/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/Melanoma_LATAM_Manuscript/Genome_Medicine_2026/Manuscript_Chakraborty_et_al_2026/Melanoma_Tirosch_cell_Annotation_CIBERSORTx/CIBERSORTX_input')

# prepare for CIBERSORTx input
write.table(
  Mixture,
  file = "Mixture.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

################################## CIBERSORTx output ##############################################################
data <- read.csv('/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/Melanoma_LATAM_Manuscript/Genome_Medicine_2026/Manuscript_Chakraborty_et_al_2026/Melanoma_Tirosch_cell_Annotation_CIBERSORTx/CIBERSORTx_output/CIBERSORTx_Job4_Results.csv')
meta <- read.xlsx('/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/Melanoma_LATAM_Manuscript/Genome_Medicine_2026/Manuscript_Chakraborty_et_al_2026/Melanoma_Tirosch_cell_Annotation_CIBERSORTx/CIBERSORTx_output/sig_mat_meta.xlsx')

#match
# match indices
m_idx <- match(meta$Gene, data$Mixture)

# columns to transfer
cols <- c("Exhausted_CD8", "Melanocytic", "Neural", "Transitory",
          "Activated_CD8", "CAF", "IFN_TAM", "LA_TAM")

# add to meta
meta[cols] <- data[m_idx, cols]
colnames(meta)[1] <- "sample_id"


#svve
write.xlsx(meta, "/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/Melanoma_LATAM_Manuscript/Genome_Medicine_2026/Manuscript_Chakraborty_et_al_2026/Melanoma_Tirosch_cell_Annotation_CIBERSORTx/CIBERSORTx_output/CIBERSORTx_output.xlsx")

#load
meta <- read.xlsx('/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/Melanoma_LATAM_Manuscript/Genome_Medicine_2026/Manuscript_Chakraborty_et_al_2026/Melanoma_Tirosch_cell_Annotation_CIBERSORTx/CIBERSORTx_output/CIBERSORTx_output.xlsx')
# selecte columns
cols <- c("Exhausted_CD8", "Melanocytic", "Neural", "Transitory",
          "Activated_CD8", "CAF", "IFN_TAM", "LA_TAM")

# ensure Response is character/factor
meta$Response <- as.character(meta$Response)

# calculate means per group
avg_NR <- colMeans(meta[meta$Response == "NR", cols], na.rm = TRUE)
avg_R  <- colMeans(meta[meta$Response == "R", cols], na.rm = TRUE)

# combine into one table
avg_df <- data.frame(
  CellType = cols,
  NR = avg_NR,
  R = avg_R
)

avg_df

########### compare with previous metadata
# first check response
# clean trailing spaces
meta$Detailed.Response <- trimws(meta$Detailed.Response)

# derive expected Response from Detailed.Response
meta$Response_expected <- ifelse(
  meta$Detailed.Response %in% c("PD", "SD"), "NR",
  ifelse(meta$Detailed.Response %in% c("CR", "PR"), "R", NA)
)

# check mismatches
mismatch <- meta[
  !is.na(meta$Response_expected) &
    meta$Response != meta$Response_expected,
]

nrow(mismatch)
mismatch

########################### Stat and visualization ############################
#load data
meta <- read.xlsx('/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/Melanoma_LATAM_Manuscript/Genome_Medicine_2026/Manuscript_Chakraborty_et_al_2026/Melanoma_Tirosch_cell_Annotation_CIBERSORTx/CIBERSORTx_output/CIBERSORTx_output.xlsx')

library(dplyr)
library(tidyr)

meta2 <- meta %>%
  mutate(
    Melanoma_total = Melanocytic + Neural + Transitory,
    CD8_total = Exhausted_CD8 + Activated_CD8,
    TAM_total = IFN_TAM + LA_TAM
  )

broad_colors <- c(
  Melanoma_total = "#E4C86E",
  CD8_total = "#67A3CC",
  TAM_total = "#FF0048",
  CAF = "#E6B0B0"
)

melanoma_colors <- c(
  Melanocytic = "#fff2cc",
  Neural = "#ffd700",
  Transitory = "#461a11"
)

tam_colors <- c(
  IFN_TAM = "#ffbedb",
  LA_TAM = "#FF0048"
)

cd8_colors <- c(
  Exhausted_CD8 = "#8bb8bc",
  Activated_CD8 = "#1850c3"
)

broad_df <- meta2 %>%
  select(Response, Melanoma_total, CD8_total, TAM_total, CAF) %>%
  pivot_longer(-Response, names_to = "CellType", values_to = "Fraction")

p1 <- ggplot(broad_df, aes(x = Response, y = Fraction, fill = CellType)) +
  stat_summary(fun = mean, geom = "bar", position = "fill") +
  scale_fill_manual(values = broad_colors) +
  scale_y_continuous(labels = scales::percent) +
  ggtitle("Broad") +
  theme_classic()

mel_df <- meta2 %>%
  select(Response, Melanocytic, Neural, Transitory) %>%
  pivot_longer(-Response, names_to = "CellType", values_to = "Fraction")

p2 <- ggplot(mel_df, aes(x = Response, y = Fraction, fill = CellType)) +
  stat_summary(fun = mean, geom = "bar", position = "fill") +
  scale_fill_manual(values = melanoma_colors) +
  ggtitle("Melanoma") +
  theme_classic()

tam_df <- meta2 %>%
  select(Response, IFN_TAM, LA_TAM) %>%
  pivot_longer(-Response, names_to = "CellType", values_to = "Fraction")

p3 <- ggplot(tam_df, aes(x = Response, y = Fraction, fill = CellType)) +
  stat_summary(fun = mean, geom = "bar", position = "fill") +
  scale_fill_manual(values = tam_colors) +
  ggtitle("TAM") +
  theme_classic()

cd8_df <- meta2 %>%
  select(Response, Exhausted_CD8, Activated_CD8) %>%
  pivot_longer(-Response, names_to = "CellType", values_to = "Fraction")

p4 <- ggplot(cd8_df, aes(x = Response, y = Fraction, fill = CellType)) +
  stat_summary(fun = mean, geom = "bar", position = "fill") +
  scale_fill_manual(values = cd8_colors) +
  ggtitle("CD8 T cells") +
  theme_classic()

library(patchwork)

final_plot <- p1 + p2 + p3 + p4 + plot_layout(ncol = 4)

final_plot

cell_types <- c(
  "Melanoma_total",
  "CD8_total",
  "TAM_total",
  "CAF",
  "Melanocytic",
  "Neural",
  "Transitory",
  "IFN_TAM",
  "LA_TAM",
  "Exhausted_CD8",
  "Activated_CD8"
)

# run tests
pvals <- sapply(cell_types, function(ct) {
  wilcox.test(meta2[[ct]] ~ meta2$Response)$p.value
})
pvals_adj <- p.adjust(pvals, method = "BH")
stat_df <- data.frame(
  CellType = cell_types,
  p_value = pvals,
  p_adj = pvals_adj
)

stat_df

write.xlsx(stat_df, '/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/Melanoma_LATAM_Manuscript/Genome_Medicine_2026/Manuscript_Chakraborty_et_al_2026/Melanoma_Tirosch_cell_Annotation_CIBERSORTx/CIBERSORTx_output/stat.xlsx')

#################################### Survival curve ####################################
# choose quantile cutoff (e.g., 75% high, 25% low)
q_high <- 0.60
q_low  <- 0.30

LA_high <- quantile(meta$LA_TAM, q_high, na.rm = TRUE)
LA_low  <- quantile(meta$LA_TAM, q_low, na.rm = TRUE)

IFN_high <- quantile(meta$IFN_TAM, q_high, na.rm = TRUE)
IFN_low  <- quantile(meta$IFN_TAM, q_low, na.rm = TRUE)

meta$TAM_group <- "Mixed-TAM"

# LA-TAM high: high LA + low IFN
meta$TAM_group[
  meta$LA_TAM >= LA_high & meta$IFN_TAM <= IFN_low
] <- "LA_TAMhigh"

# IFN-TAM high: high IFN + low LA
meta$TAM_group[
  meta$IFN_TAM >= IFN_high & meta$LA_TAM <= LA_low
] <- "IFN_TAMhigh"

table(meta$TAM_group)

library(survival)
library(survminer)

survInput <- meta[, c(
  "sample_id", "Overall_survival", "Alive", "Therapy", "Cohort", "Response", "TAM_group"
)]

# time
survInput$TIME <- as.numeric(survInput$Overall_survival)

# status: 1 = event/death, 0 = censored/alive
survInput$STATUS <- ifelse(
  survInput$Alive %in% c("NO", "No", "no", 0, "0"),
  1, 0
)

# drop missing values
survInput <- survInput[!is.na(survInput$TIME) & !is.na(survInput$STATUS) & !is.na(survInput$TAM_group), ]

# optional: set factor order
survInput$TAM_group <- factor(
  survInput$TAM_group,
  levels = c("IFN_TAMhigh", "Mixed-TAM", "LA_TAMhigh")
)

fit <- survfit(Surv(TIME, STATUS) ~ TAM_group, data = survInput)

p <- ggsurvplot(
  fit,
  data = survInput,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  xlim = c(0, 2000),
  xlab = "Time in days",
  break.time.by = 500,
  ggtheme = theme_light(),
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE
)

p

survInput$TAM_group <- factor(
  survInput$TAM_group,
  levels = c("IFN_TAMhigh", "Mixed-TAM", "LA_TAMhigh")
)

fit <- survfit(Surv(TIME, STATUS) ~ TAM_group, data = survInput)

tam_colors <- c("#ffbedb", "#E8DBB3", "#FF0048")

p <- ggsurvplot(
  fit,
  data = survInput,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  xlim = c(0, 2000),
  xlab = "Time in days",
  break.time.by = 500,
  palette = tam_colors,
  ggtheme = theme_light(),
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE
)

p

########################### COXPH Uni and Multivariate analysis ######################
########################### COXPH Uni and Multivariate analysis ##########################

library(dplyr)
library(survival)
library(survminer)
library(openxlsx)

# Load meta
meta <- read.xlsx("/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/Melanoma_LATAM_Manuscript/Genome_Medicine_2026/Manuscript_Chakraborty_et_al_2026/Melanoma_Tirosch_cell_Annotation_CIBERSORTx/CIBERSORTx_output/CIBERSORTx_output.xlsx")

# Variables of interest
cell_types <- c(
  "Exhausted_CD8", "Activated_CD8",
  "Melanocytic", "Neural", "Transitory",
  "CAF", "IFN_TAM", "LA_TAM"
)

# Build survival input
survInput <- meta[, c(
  "sample_id", "Overall_survival", "Alive", "Cohort", "Therapy", "Response",
  cell_types
)]

# Time and status
survInput$TIME <- as.numeric(survInput$Overall_survival)

survInput$STATUS <- ifelse(
  survInput$Alive %in% c("NO", "No", "no", 0, "0"),
  1, 0
)

# Factors
survInput$Response <- factor(survInput$Response, levels = c("R", "NR"))
survInput$Therapy <- factor(survInput$Therapy)
survInput$Cohort <- factor(survInput$Cohort)

# Remove missing time/status and missing predictor rows
survInput <- survInput %>%
  filter(!is.na(TIME), !is.na(STATUS))

# Optional but recommended:
# scale continuous predictors so HR is per 1 SD increase
survInput_scaled <- survInput
survInput_scaled[cell_types] <- scale(survInput_scaled[cell_types])

# ----------------------------- UNIVARIATE COX ------------------------------------------

univ_results <- lapply(cell_types, function(var) {
  
  df_sub <- survInput_scaled %>%
    select(TIME, STATUS, all_of(var)) %>%
    filter(complete.cases(.))
  
  form <- as.formula(paste0("Surv(TIME, STATUS) ~ ", var))
  fit <- coxph(form, data = df_sub)
  sm <- summary(fit)
  
  ph <- cox.zph(fit)
  ph_p <- ph$table[1, "p"]
  
  data.frame(
    CellType = var,
    N = nrow(df_sub),
    Events = sum(df_sub$STATUS, na.rm = TRUE),
    beta = as.numeric(sm$coefficients[1, "coef"]),
    HR = as.numeric(sm$coefficients[1, "exp(coef)"]),
    lower95 = as.numeric(sm$conf.int[1, "lower .95"]),
    upper95 = as.numeric(sm$conf.int[1, "upper .95"]),
    p.value = as.numeric(sm$coefficients[1, "Pr(>|z|)"]),
    PH_p.value = as.numeric(ph_p)
  )
})

univ_res <- bind_rows(univ_results)

# Multiple testing correction
univ_res$FDR <- p.adjust(univ_res$p.value, method = "BH")

# Formatted HR string
univ_res$`HR (95% CI)` <- paste0(
  sprintf("%.2f", univ_res$HR), " (",
  sprintf("%.2f", univ_res$lower95), "-",
  sprintf("%.2f", univ_res$upper95), ")"
)

# Order results
univ_res <- univ_res %>%
  arrange(p.value)

# Clean display table
univ_res_table <- univ_res %>%
  select(
    CellType, N, Events, beta, HR, lower95, upper95,
    `HR (95% CI)`, p.value, FDR, PH_p.value
  )

print(univ_res_table)

# Export
setwd('/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/Melanoma_LATAM_Manuscript/Genome_Medicine_2026/Manuscript_Chakraborty_et_al_2026/Melanoma_Tirosch_cell_Annotation_CIBERSORTx/CIBERSORTx_output/Survival')

write.xlsx(
  list(
    Univariate_Cox = univ_res_table
  ),
  file = "coxph_univariate_celltypes.xlsx",
  overwrite = TRUE
)

# ----------------------------- OPTIONAL FOREST PLOT ------------------------------------

library(ggplot2)

p_univ_forest <- ggplot(
  univ_res,
  aes(x = reorder(CellType, HR), y = HR)
) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower95, ymax = upper95), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  coord_flip() +
  theme_minimal(base_size = 13) +
  labs(
    title = "Univariate Cox regression",
    x = "",
    y = "Hazard ratio per 1 SD increase"
  )

print(p_univ_forest)


#################################### MULTIVARIATE ANALYSIS ######################################
########################### COXPH Multivariate analysis #########################################

library(dplyr)
library(survival)
library(survminer)
library(openxlsx)
library(ggplot2)

# Load meta
meta <- read.xlsx("/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/Melanoma_LATAM_Manuscript/Genome_Medicine_2026/Manuscript_Chakraborty_et_al_2026/Melanoma_Tirosch_cell_Annotation_CIBERSORTx/CIBERSORTx_output/CIBERSORTx_output.xlsx")

# Variables of interest
cell_types <- c(
  "Transitory",
  "IFN_TAM",
  "Activated_CD8"
)

# Build survival input
survInput <- meta[, c(
  "sample_id", "Overall_survival", "Alive", "Cohort", "Therapy", "Response",
  cell_types
)]

# Time and status
survInput$TIME <- as.numeric(survInput$Overall_survival)

survInput$STATUS <- ifelse(
  survInput$Alive %in% c("NO", "No", "no", 0, "0"),
  1, 0
)

# Factors
survInput$Response <- factor(survInput$Response, levels = c("R", "NR"))
survInput$Therapy <- factor(survInput$Therapy)
survInput$Cohort <- factor(survInput$Cohort)

# Remove missing values
survInput <- survInput %>%
  filter(!is.na(TIME), !is.na(STATUS)) %>%
  filter(complete.cases(across(all_of(c(cell_types, "Therapy", "Cohort")))))

# Scale continuous predictors so HR is per 1 SD increase
survInput_scaled <- survInput
survInput_scaled[cell_types] <- scale(survInput_scaled[cell_types])

# ----------------------------- MULTIVARIATE COX -----------------------------------------

multi_formula <- as.formula(
  paste(
    "Surv(TIME, STATUS) ~",
    paste(c(cell_types, "Therapy", "Cohort"), collapse = " + ")
  )
)

multi_model <- coxph(multi_formula, data = survInput_scaled)
multi_summary <- summary(multi_model)

# PH test
multi_ph <- cox.zph(multi_model)
multi_ph_table <- as.data.frame(multi_ph$table)
multi_ph_table$Term <- rownames(multi_ph_table)
rownames(multi_ph_table) <- NULL

# Extract coefficient table
coef_df <- as.data.frame(multi_summary$coefficients)
conf_df <- as.data.frame(multi_summary$conf.int)

multi_res <- data.frame(
  Term = rownames(coef_df),
  beta = coef_df$coef,
  HR = conf_df$`exp(coef)`,
  lower95 = conf_df$`lower .95`,
  upper95 = conf_df$`upper .95`,
  p.value = coef_df$`Pr(>|z|)`,
  stringsAsFactors = FALSE
)

# Add PH p-values
multi_res <- multi_res %>%
  left_join(
    multi_ph_table %>% select(Term, PH_p.value = p),
    by = "Term"
  )

# Add formatted HR string
multi_res$`HR (95% CI)` <- paste0(
  sprintf("%.2f", multi_res$HR), " (",
  sprintf("%.2f", multi_res$lower95), "-",
  sprintf("%.2f", multi_res$upper95), ")"
)

# Add model-level info
multi_res$N <- multi_summary$n
multi_res$Events <- multi_summary$nevent
multi_res$FDR <- p.adjust(multi_res$p.value, method = "BH")

# Reorder columns
multi_res_table <- multi_res %>%
  filter(Term %in% cell_types) %>%   # 👈 THIS is the key line
  select(
    Term, N, Events, beta, HR, lower95, upper95,
    `HR (95% CI)`, p.value, FDR
  ) %>%
  arrange(p.value)

print(multi_res_table)

# ----------------------------- EXPORT ---------------------------------------------------

setwd("/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/Melanoma_LATAM_Manuscript/Genome_Medicine_2026/Manuscript_Chakraborty_et_al_2026/Melanoma_Tirosch_cell_Annotation_CIBERSORTx/CIBERSORTx_output")

write.xlsx(
  list(
    Multivariate_Cox = multi_res_table,
    PH_Test = multi_ph_table
  ),
  file = "coxph_multivariate_celltypes.xlsx",
  overwrite = TRUE
)

# ----------------------------- OPTIONAL FOREST PLOT ------------------------------------

# Keep only main biological variables for forest plot
plot_df <- multi_res %>%
  filter(Term %in% cell_types)

p_multi_forest <- ggplot(
  plot_df,
  aes(x = reorder(Term, HR), y = HR)
) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower95, ymax = upper95), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  coord_flip() +
  theme_minimal(base_size = 13) +
  labs(
    title = "Multivariate Cox regression",
    x = "",
    y = "Hazard ratio per 1 SD increase"
  )

print(p_multi_forest)



#########################################################################################
################# XgBoost/ RAndom Forest ################################################
#########################################################################################
#########################################################################################
################# XGBoost LOCO with 20 repeated iterations ##############################
#########################################################################################
library(dplyr)
library(tidyr)
library(xgboost)
library(pROC)
library(openxlsx)
library(ggplot2)

data <- read.xlsx('/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/Melanoma_LATAM_Manuscript/Genome_Medicine_2026/OLD2/Melanoma_Tirosch_cell_Annotation_CIBERSORTx/CIBERSORTx_output/CIBERSORTx_output.xlsx')

features <- c(
  "Exhausted_CD8",
  "Melanocytic",
  "Neural",
  "Transitory",
  "Activated_CD8",
  "CAF",
  "IFN_TAM",
  "LA_TAM"
)

data2 <- data %>%
  mutate(Response_bin = factor(Response, levels = c("NR", "R"))) %>%
  drop_na(all_of(features), Response_bin)

cohorts <- unique(data2$Cohort)

cohort_colors <- c(
  "Gide" = "#2A9D8F",
  "Hugo" = "#A8D5BA",
  "Liu" = "#E9D8A6",
  "MGH" = "#FFC300",
  "Riaz" = "#F8961E",
  "Van Allen" = "#D62828"
)

make_val_indices <- function(y, prop = 0.2, seed = 123) {
  set.seed(seed)
  y <- as.factor(y)
  idx_by_class <- split(seq_along(y), y)
  
  val_idx <- unlist(lapply(idx_by_class, function(idx) {
    n_take <- max(1, round(length(idx) * prop))
    sample(idx, size = n_take)
  }))
  
  sort(val_idx)
}

# ----------------------------- iteration settings --------------------------------------
n_iterations <- 10
iteration_seeds <- 1000 + seq_len(n_iterations)

# ----------------------------- storage objects -----------------------------------------
auc_results <- list()
best_params_all <- list()
importance_all <- list()
preds_all <- list()

# ============================= OUTER ITERATION LOOP ====================================
for (iter in seq_len(n_iterations)) {
  
  current_seed <- iteration_seeds[iter]
  
  cat("\n############################################################\n")
  cat("Iteration:", iter, "| Seed:", current_seed, "\n")
  cat("############################################################\n")
  
  for (cohort in cohorts) {
    
    cat("\n====================================================\n")
    cat("Iteration:", iter, "| Testing on cohort:", cohort, "\n")
    cat("====================================================\n")
    
    train_data_full <- data2 %>% filter(Cohort != cohort)
    test_data       <- data2 %>% filter(Cohort == cohort)
    
    x_train_full <- as.matrix(train_data_full[, features, drop = FALSE])
    y_train_full <- train_data_full$Response_bin
    
    x_test <- as.matrix(test_data[, features, drop = FALSE])
    y_test <- test_data$Response_bin
    y_test_num <- ifelse(y_test == "R", 1, 0)
    
    val_idx <- make_val_indices(y_train_full, prop = 0.4, seed = current_seed)
    train_idx <- setdiff(seq_len(nrow(train_data_full)), val_idx)
    
    x_inner_train <- x_train_full[train_idx, , drop = FALSE]
    y_inner_train <- y_train_full[train_idx]
    
    x_val <- x_train_full[val_idx, , drop = FALSE]
    y_val <- y_train_full[val_idx]
    y_val_num <- ifelse(y_val == "R", 1, 0)
    
    tune_grid <- expand.grid(
      nrounds = c(50, 100, 200, 300),
      learning_rate = c(0.01, 0.05, 0.1),
      max_depth = c(2, 3, 4),
      min_child_weight = c(1, 3),
      subsample = c(0.8, 1.0),
      colsample_bytree = c(0.8, 1.0),
      gamma = c(0, 1)
    )
    
    tuning_results <- list()
    counter <- 1
    
    for (i in seq_len(nrow(tune_grid))) {
      
      fit <- tryCatch(
        {
          xgboost(
            x = x_inner_train,
            y = y_inner_train,
            objective = "binary:logistic",
            eval_metric = "auc",
            nrounds = tune_grid$nrounds[i],
            max_depth = tune_grid$max_depth[i],
            learning_rate = tune_grid$learning_rate[i],
            min_child_weight = tune_grid$min_child_weight[i],
            subsample = tune_grid$subsample[i],
            colsample_bytree = tune_grid$colsample_bytree[i],
            gamma = tune_grid$gamma[i],
            verbosity = 0,
            seed = current_seed
          )
        },
        error = function(e) {
          message("Tuning error in iteration ", iter, ", cohort ", cohort, ", row ", i, ": ", e$message)
          NULL
        }
      )
      
      if (is.null(fit)) next
      
      preds_val <- tryCatch(
        {
          predict(fit, x_val, type = "response")
        },
        error = function(e) NULL
      )
      
      if (is.null(preds_val)) next
      if (length(unique(y_val_num)) < 2) next
      
      val_auc <- tryCatch(
        {
          as.numeric(auc(roc(response = y_val_num, predictor = preds_val, quiet = TRUE)))
        },
        error = function(e) NA_real_
      )
      
      if (!is.finite(val_auc)) next
      
      tuning_results[[counter]] <- data.frame(
        Iteration = iter,
        Cohort = cohort,
        nrounds = tune_grid$nrounds[i],
        learning_rate = tune_grid$learning_rate[i],
        max_depth = tune_grid$max_depth[i],
        min_child_weight = tune_grid$min_child_weight[i],
        subsample = tune_grid$subsample[i],
        colsample_bytree = tune_grid$colsample_bytree[i],
        gamma = tune_grid$gamma[i],
        val_auc = val_auc
      )
      
      counter <- counter + 1
    }
    
    tuning_results_df <- bind_rows(tuning_results)
    
    if (nrow(tuning_results_df) == 0) {
      cat("No valid tuning results for", cohort, "in iteration", iter, "- using default parameters.\n")
      best_row <- data.frame(
        Iteration = iter,
        Cohort = cohort,
        nrounds = 200,
        learning_rate = 0.05,
        max_depth = 3,
        min_child_weight = 1,
        subsample = 0.8,
        colsample_bytree = 0.8,
        gamma = 0,
        val_auc = NA_real_
      )
    } else {
      best_row <- tuning_results_df %>%
        arrange(desc(val_auc)) %>%
        slice(1)
      cat("Best validation AUC:", round(best_row$val_auc, 3), "\n")
    }
    
    best_params_all[[length(best_params_all) + 1]] <- best_row
    
    final_model <- xgboost(
      x = x_train_full,
      y = y_train_full,
      objective = "binary:logistic",
      eval_metric = "auc",
      nrounds = best_row$nrounds,
      max_depth = best_row$max_depth,
      learning_rate = best_row$learning_rate,
      min_child_weight = best_row$min_child_weight,
      subsample = best_row$subsample,
      colsample_bytree = best_row$colsample_bytree,
      gamma = best_row$gamma,
      verbosity = 0,
      seed = current_seed
    )
    
    preds <- predict(final_model, x_test, type = "response")
    
    roc_obj <- roc(response = y_test_num, predictor = preds, quiet = TRUE)
    auc_val <- as.numeric(auc(roc_obj))
    
    cat("Test AUC:", round(auc_val, 3), "\n")
    
    auc_results[[length(auc_results) + 1]] <- data.frame(
      Iteration = iter,
      Cohort = cohort,
      AUC = auc_val,
      n_test = nrow(test_data)
    )
    
    preds_all[[length(preds_all) + 1]] <- data.frame(
      Iteration = iter,
      Cohort = cohort,
      Response = y_test_num,
      Prediction = preds
    )
    
    imp <- xgb.importance(feature_names = features, model = final_model)
    imp$Iteration <- iter
    imp$Cohort <- cohort
    importance_all[[length(importance_all) + 1]] <- imp
  }
}

# ----------------------------- combine results -----------------------------------------
auc_results_df <- bind_rows(auc_results)
best_params_df <- bind_rows(best_params_all)
importance_df <- bind_rows(importance_all)
all_preds_df <- bind_rows(preds_all)

# ----------------------------- summary: per cohort -------------------------------------
auc_summary_cohort <- auc_results_df %>%
  group_by(Cohort) %>%
  summarise(
    n_iterations = n(),
    Mean_AUC = mean(AUC, na.rm = TRUE),
    SD_AUC = sd(AUC, na.rm = TRUE),
    SEM_AUC = SD_AUC / sqrt(n_iterations),
    Min_AUC = min(AUC, na.rm = TRUE),
    Max_AUC = max(AUC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(match(Cohort, cohorts))

print(auc_summary_cohort)

# ----------------------------- summary: overall across fold x iteration ----------------
auc_summary_overall <- auc_results_df %>%
  summarise(
    n_total = n(),
    Mean_AUC = mean(AUC, na.rm = TRUE),
    SD_AUC = sd(AUC, na.rm = TRUE),
    SEM_AUC = SD_AUC / sqrt(n_total),
    Min_AUC = min(AUC, na.rm = TRUE),
    Max_AUC = max(AUC, na.rm = TRUE)
  )

print(auc_summary_overall)

# ----------------------------- optional: overall mean by iteration ---------------------
auc_summary_iteration <- auc_results_df %>%
  group_by(Iteration) %>%
  summarise(
    Mean_AUC = mean(AUC, na.rm = TRUE),
    SD_AUC = sd(AUC, na.rm = TRUE),
    .groups = "drop"
  )

print(auc_summary_iteration)

# ----------------------------- ROC summary plot ----------------------------------------
p_auc_iter <- ggplot(auc_results_df, aes(x = Cohort, y = AUC, color = Cohort)) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 2) +
  stat_summary(fun = mean, geom = "point", size = 4, shape = 18, color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "black") +
  scale_color_manual(values = cohort_colors) +
  theme_minimal(base_size = 13) +
  labs(
    title = "AUC across 10 repeated iterations",
    x = "Test cohort",
    y = "AUC"
  )

print(p_auc_iter)

# ----------------------------- feature importance summary ------------------------------
importance_summary <- importance_df %>%
  group_by(Feature) %>%
  summarise(
    MeanGain = mean(Gain, na.rm = TRUE),
    SDGain = sd(Gain, na.rm = TRUE),
    SEMGain = SDGain / sqrt(n()),
    .groups = "drop"
  ) %>%
  arrange(desc(MeanGain))

print(importance_summary)

p_imp <- ggplot(
  importance_summary,
  aes(x = reorder(Feature, MeanGain), y = MeanGain)
) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(aes(ymin = MeanGain - SEMGain, ymax = MeanGain + SEMGain), width = 0.2) +
  coord_flip() +
  theme_minimal(base_size = 13) +
  ylab("Mean Importance (Gain) ± SEM") +
  xlab("") +
  ggtitle("Feature Importance Across Cohorts and Iterations")

print(p_imp)

# ----------------------------- pooled confusion-ready predictions ----------------------
all_preds_df <- all_preds_df %>%
  mutate(
    Predicted_class = ifelse(Prediction >= 0.5, 1, 0),
    Outcome = case_when(
      Response == 1 & Predicted_class == 1 ~ "TP",
      Response == 0 & Predicted_class == 0 ~ "TN",
      Response == 0 & Predicted_class == 1 ~ "FP",
      Response == 1 & Predicted_class == 0 ~ "FN"
    )
  )

setwd('/Users/sajibchakraborty/Documents/EkoEstMed/Melanoma Project/Melanoma_LATAM_Manuscript/Genome_Medicine_2026/Manuscript_Chakraborty_et_al_2026/Melanoma_Tirosch_cell_Annotation_CIBERSORTx/CIBERSORTx_output/XGboost/HyParam_Tuning')

# ----------------------------- save if needed ------------------------------------------
write.csv(auc_results_df, "AUC_results_10_iterations.csv", row.names = FALSE)
write.csv(auc_summary_cohort, "AUC_summary_by_cohort_10_iterations.csv", row.names = FALSE)
write.csv(importance_summary, "Importance_summary_10_iterations.csv", row.names = FALSE)
write.csv(best_params_df, "Best_params_10_iterations.csv", row.names = FALSE)
write.csv(all_preds_df, "Predictions_10_iterations.csv", row.names = FALSE)

library(openxlsx)

write.xlsx(
  list(
    AUC_results = auc_results_df,
    AUC_summary = auc_summary_cohort,
    Importance = importance_summary,
    Best_params = best_params_df,
    Predictions = all_preds_df
  ),
  file = "XGBoost_LOCO_10_iterations_results.xlsx",
  overwrite = TRUE
)

# ROC plot 
roc_list <- list()
auc_summary <- data.frame()

for (cohort in unique(all_preds_df$Cohort)) {
  
  df_sub <- all_preds_df %>% filter(Cohort == cohort)
  
  roc_obj <- roc(
    response = df_sub$Response,
    predictor = df_sub$Prediction,
    quiet = TRUE
  )
  
  roc_list[[cohort]] <- roc_obj
  
  # mean AUC across iterations
  auc_iter <- auc_results_df %>% 
    filter(Cohort == cohort) %>% 
    pull(AUC)
  
  auc_summary <- rbind(
    auc_summary,
    data.frame(
      Cohort = cohort,
      Mean_AUC = mean(auc_iter),
      SEM = sd(auc_iter)/sqrt(length(auc_iter))
    )
  )
}

# plot
roc_df_list <- list()

for (cohort in names(roc_list)) {
  roc_df_list[[cohort]] <- data.frame(
    Specificity = roc_list[[cohort]]$specificities,
    Sensitivity = roc_list[[cohort]]$sensitivities,
    Cohort = cohort
  )
}

roc_df <- bind_rows(roc_df_list)

# labels
cohort_labels <- setNames(
  paste0(
    auc_summary$Cohort,
    " (AUC = ",
    round(auc_summary$Mean_AUC, 3),
    " ± ",
    round(auc_summary$SEM, 3),
    ")"
  ),
  auc_summary$Cohort
)

ggplot(roc_df, aes(x = 1 - Specificity, y = Sensitivity, color = Cohort)) +
  geom_line(linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = cohort_colors,
                     labels = cohort_labels) +
  theme_minimal(base_size = 13) +
  labs(
    title = "ROC curves (mean across 10 iterations)",
    x = "False Positive Rate",
    y = "True Positive Rate",
    color = "Test cohort"
  )

############################# Contengency plot of prediction test cohort #################
############################# Contingency plot of prediction by test cohort #############
library(dplyr)
library(tidyr)
library(ggplot2)
library(colorspace)

# Base cohort colors
cohort_colors <- c(
  "Gide" = "#2A9D8F",
  "Hugo" = "#A8D5BA",
  "Liu" = "#E9D8A6",
  "MGH" = "#FFC300",
  "Riaz" = "#F8961E",
  "Van Allen" = "#D62828"
)

# Use repeated-iteration predictions table from the 10-iteration pipeline:
# all_preds_df columns:
# Iteration | Cohort | Response | Prediction

all_preds2 <- all_preds_df %>%
  mutate(
    Predicted_class = ifelse(Prediction >= 0.5, 1, 0),
    Actual = ifelse(Response == 1, "R", "NR"),
    Predicted = ifelse(Predicted_class == 1, "R", "NR"),
    Type = case_when(
      Actual == "R"  & Predicted == "R"  ~ "TP",
      Actual == "NR" & Predicted == "NR" ~ "TN",
      Actual == "NR" & Predicted == "R"  ~ "FP",
      Actual == "R"  & Predicted == "NR" ~ "FN"
    )
  )

# Count TP/TN/FP/FN within each iteration and cohort
conf_iter <- all_preds2 %>%
  count(Iteration, Cohort, Actual, Predicted, Type, .drop = FALSE)

# Ensure all 4 cells exist for every iteration x cohort
conf_iter <- conf_iter %>%
  complete(
    Iteration,
    Cohort,
    Actual = c("NR", "R"),
    Predicted = c("NR", "R"),
    fill = list(n = 0)
  ) %>%
  mutate(
    Type = case_when(
      Actual == "R"  & Predicted == "R"  ~ "TP",
      Actual == "NR" & Predicted == "NR" ~ "TN",
      Actual == "NR" & Predicted == "R"  ~ "FP",
      Actual == "R"  & Predicted == "NR" ~ "FN"
    )
  )

# Add per-iteration percentages within each cohort
conf_iter <- conf_iter %>%
  group_by(Iteration, Cohort) %>%
  mutate(
    Total = sum(n),
    Percent = 100 * n / Total
  ) %>%
  ungroup()

# Summarise across iterations
conf_summary <- conf_iter %>%
  group_by(Cohort, Actual, Predicted, Type) %>%
  summarise(
    Mean_n = mean(n),
    SD_n = sd(n),
    SEM_n = SD_n / sqrt(n()),
    Mean_percent = mean(Percent),
    SD_percent = sd(Percent),
    SEM_percent = SD_percent / sqrt(n()),
    .groups = "drop"
  )

# Make shades for each cohort
make_shades <- function(base_col) {
  c(
    "TP" = darken(base_col, 0.15),
    "TN" = lighten(base_col, 0.10),
    "FP" = lighten(base_col, 0.35),
    "FN" = darken(base_col, 0.35)
  )
}

shade_df <- do.call(
  rbind,
  lapply(names(cohort_colors), function(coh) {
    s <- make_shades(cohort_colors[[coh]])
    data.frame(
      Cohort = coh,
      Type = names(s),
      Fill = unname(s),
      stringsAsFactors = FALSE
    )
  })
)

# Join colors and labels
conf_plot_df <- conf_summary %>%
  left_join(shade_df, by = c("Cohort", "Type")) %>%
  mutate(
    Actual = factor(Actual, levels = c("NR", "R")),
    Predicted = factor(Predicted, levels = c("NR", "R")),
    Label = paste0(
      Type,
      "\nMean=", round(Mean_n, 1),
      "\nSEM=", round(SEM_n, 1),
      "\n", round(Mean_percent, 1), "%"
    )
  )

# Plot faceted 2x2 confusion matrices
ggplot(conf_plot_df, aes(x = Predicted, y = Actual)) +
  geom_tile(aes(fill = Fill), color = "white", linewidth = 1.2) +
  geom_text(aes(label = Label), fontface = "bold", size = 3.8, lineheight = 0.95) +
  scale_fill_identity() +
  facet_wrap(~ Cohort, ncol = 3) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Confusion matrices by test cohort\n(mean across iterations)",
    x = "Predicted",
    y = "Actual"
  ) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text = element_text(face = "bold")
  )

############################# Calculate neighborhood fractions ###########################
library(dplyr)
library(tidyr)
library(ggplot2)

library(dplyr)

df <- data %>%
  filter(Detailed.Response %in% c("CR", "PR", "SD", "PD")) %>%
  mutate(
    LA_TAM_neighborhood = LA_TAM + Neural + Exhausted_CD8,
    IFN_TAM_neighborhood = IFN_TAM + Exhausted_CD8
  ) %>%
  group_by(sample_id, Detailed.Response) %>%
  summarise(
    LA = mean(LA_TAM_neighborhood, na.rm = TRUE),
    IFN = mean(IFN_TAM_neighborhood, na.rm = TRUE),
    .groups = "drop"
  )

# order
df$Detailed.Response <- factor(df$Detailed.Response,
                               levels = c("CR", "PR", "SD", "PD"))

summary_df <- df %>%
  group_by(Detailed.Response) %>%
  summarise(
    n = n(),
    LA_mean = -mean(LA),   # 🔥 NEGATIVE for left side
    LA_sem = sd(LA) / sqrt(n),
    IFN_mean = mean(IFN),
    IFN_sem = sd(IFN) / sqrt(n)
  )

library(ggplot2)

ggplot(summary_df) +
  
  # LA-TAM (left)
  geom_bar(aes(x = Detailed.Response, y = LA_mean),
           stat = "identity",
           fill = "#FF2D1A",
           width = 0.6) +
  
  geom_errorbar(aes(x = Detailed.Response,
                    ymin = LA_mean - LA_sem,
                    ymax = LA_mean + LA_sem),
                width = 0.2) +
  
  # IFN-TAM (right)
  geom_bar(aes(x = Detailed.Response, y = IFN_mean),
           stat = "identity",
           fill = "#5B3FD3",
           width = 0.6) +
  
  geom_errorbar(aes(x = Detailed.Response,
                    ymin = IFN_mean - IFN_sem,
                    ymax = IFN_mean + IFN_sem),
                width = 0.2) +
  
  coord_flip() +
  
  geom_vline(xintercept = 0, color = "black") +
  
  theme_minimal(base_size = 14) +
  labs(
    x = "",
    y = "Neighborhood fraction",
    title = "TAM neighborhood continuum"
  )

df %>%
  count(Detailed.Response) %>%
  mutate(label = paste0(Detailed.Response, " (n=", n, ")"))

################################### Final ROC curve #########################################
library(dplyr)
library(xgboost)
library(pROC)
library(ggplot2)
library(openxlsx)

# Load data
data <- read.xlsx("CIBERSORTx_output.xlsx")

features <- c(
  "Exhausted_CD8", "Melanocytic", "Neural", "Transitory",
  "Activated_CD8", "CAF", "IFN_TAM", "LA_TAM"
)

data2 <- data %>%
  mutate(Response_bin = factor(Response, levels = c("NR", "R"))) %>%
  drop_na(all_of(features), Response_bin)

cohorts <- unique(data2$Cohort)

roc_list <- list()
auc_results <- data.frame()

# ===================== LOCO (single run) =====================
for (cohort in cohorts) {
  
  cat("Testing on cohort:", cohort, "\n")
  
  train_data <- data2 %>% filter(Cohort != cohort)
  test_data  <- data2 %>% filter(Cohort == cohort)
  
  x_train <- as.matrix(train_data[, features])
  y_train <- train_data$Response_bin
  
  x_test <- as.matrix(test_data[, features])
  y_test <- test_data$Response_bin
  y_test_num <- ifelse(y_test == "R", 1, 0)
  
  # ---- model (fixed parameters, no tuning) ----
  model <- xgboost(
    x = x_train,
    y = y_train,
    objective = "binary:logistic",
    eval_metric = "auc",
    nrounds = 200,
    max_depth = 3,
    learning_rate = 0.05,
    subsample = 0.8,
    colsample_bytree = 0.8,
    verbosity = 0,
    seed = 123
  )
  
  preds <- predict(model, x_test)
  
  roc_obj <- roc(response = y_test_num, predictor = preds, quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  
  cat("AUC:", round(auc_val, 3), "\n")
  
  roc_list[[cohort]] <- roc_obj
  
  auc_results <- rbind(
    auc_results,
    data.frame(Cohort = cohort, AUC = auc_val)
  )
}

# ===================== Build ROC dataframe =====================
roc_df_list <- list()

for (cohort in names(roc_list)) {
  roc_df_list[[cohort]] <- data.frame(
    Specificity = roc_list[[cohort]]$specificities,
    Sensitivity = roc_list[[cohort]]$sensitivities,
    Cohort = cohort
  )
}

roc_df <- bind_rows(roc_df_list)

# Labels with AUC
cohort_labels <- setNames(
  paste0(auc_results$Cohort, " (AUC = ", round(auc_results$AUC, 3), ")"),
  auc_results$Cohort
)

# Colors
cohort_colors <- c(
  "Gide" = "#2A9D8F",
  "Hugo" = "#A8D5BA",
  "Liu" = "#E9D8A6",
  "MGH" = "#FFC300",
  "Riaz" = "#F8961E",
  "Van Allen" = "#D62828"
)

# ===================== Plot =====================
ggplot(roc_df, aes(x = 1 - Specificity, y = Sensitivity, color = Cohort)) +
  geom_line(linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = cohort_colors, labels = cohort_labels) +
  theme_minimal(base_size = 13) +
  labs(
    title = "ROC curves (single LOCO run)",
    x = "False Positive Rate",
    y = "True Positive Rate",
    color = "Test cohort"
  )

# ===================== Plot2 =====================
roc_df_list <- list()

for (cohort in names(roc_list)) {
  
  roc_df_list[[cohort]] <- data.frame(
    FPR = 1 - roc_list[[cohort]]$specificities,
    TPR = roc_list[[cohort]]$sensitivities,
    Cohort = cohort
  ) %>%
    dplyr::arrange(FPR, TPR)
}

roc_df <- dplyr::bind_rows(roc_df_list)

ggplot(roc_df, aes(x = FPR, y = TPR, color = Cohort)) +
  geom_step(linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = cohort_colors,
                     labels = cohort_labels) +
  theme_minimal(base_size = 13) +
  labs(
    title = "ROC curves",
    x = "False Positive Rate",
    y = "True Positive Rate",
    color = "Test cohort"
  )
