# Libraries ----
library(tidyverse)
library(sleuth)
library(biomaRt)
library(gumshoe)
library(ggrepel)
library(multcompView)
library(rcompanion)
library(ggpubr)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
setwd("~/tutorials/Tutorial_Pseudo_Bulk")

# Metadata File ----
## Load Metadata File ----
metadata <- read.table("metadata.txt", header = TRUE)
## Re-level Metadata ----
# Convert the tissue type column to a factor and, based on how R works, NF will be the first factor level.
metadata$tissue <- as.factor(metadata$tissue)

# Create a new variable that will contain a relevled version of the metadata
metadata_releveled <- metadata
metadata_releveled$tissue <- relevel(metadata_releveled$tissue, ref = "TH")

# Gumshoe Data Frame Creation ----
## Analysis Information ----
# Define the metadata file name(s), associated model name(s), and the corresponding model(s) to be used.
# NOTE: The order of the model data must be align with the model names.
metadata_names <- c("metadata_nominal_NF_baseline",
                    "metadata_nominal_TH_baseline")
model_names <- c("NF_model_sex,NF_model_sex_tissue,NF_model_interaction",
                 "NF_t2g_model_sex,NF_t2g_model_sex_tissue,NF_t2g_model_interaction",
                 "TH_model_sex,TH_model_sex_tissue,TH_model_interaction",
                 "TH_t2g_model_sex,TH_t2g_model_sex_tissue,TH_t2g_model_interaction")
model_data <- c("~sex, ~sex + tissue, ~sex*tissue")
model_parameters <- c("",
                      "target_map = t2g, aggregation_column = 'ext_gene', gene_mode = TRUE")
## Data Frame Creation ----
# Take the metadata names, model names, and model data information and combine it into a single data frame.
analysis_data <- data.frame(metadata_name = metadata_names, 
                            metadata_file = tibble(list(metadata)), 
                            model_name = model_names, 
                            model_data = model_data, 
                            model_parameters = model_parameters)

# Replace the metadata file with the re-leveled version.
analysis_data[[2]][[2]] <- metadata_releveled
analysis_data[[2]][[4]] <- metadata_releveled

levels(analysis_data[[2]][[1]]$tissue)
# "NF"  "NP"  "PEP" "TH"
levels(analysis_data[[2]][[2]]$tissue)
# "TH"  "NF"  "NP"  "PEP"
levels(analysis_data[[2]][[3]]$tissue)
# "NF"  "NP"  "PEP" "TH"
levels(analysis_data[[2]][[4]]$tissue)
# "TH"  "NF"  "NP"  "PEP"

# Gumshoe Function Calls ----
# To get an association between the transcript versions and the gene ids 
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
t2g <- getBM(attributes = c("ensembl_transcript_id_version", "ensembl_gene_id","external_gene_name"), mart = mart)
colnames(t2g) <- c("target_id", "ens_gene", "ext_gene")

sleuth_interpret(analysis_data, num_core = 12)

# Check the model RSS values ----
rss_df <- data.frame(models = c("so_NF_model_sex", "so_NF_t2g_model_sex", "so_TH_model_sex", "so_TH_t2g_model_sex",
                                "so_NF_model_sex_tissue", "so_NF_t2g_model_sex_tissue", "so_TH_model_sex_tissue", "so_TH_t2g_model_sex_tissue",
                                "so_NF_model_interaction", "so_TH_model_interaction", "so_NF_t2g_model_interaction", "so_TH_t2g_model_interaction"),
                     model_rss = c(sleuth_model_rss(so_NF_model_sex), sleuth_model_rss(so_TH_model_sex), sleuth_model_rss(so_NF_t2g_model_sex),  sleuth_model_rss(so_TH_t2g_model_sex),
                                   sleuth_model_rss(so_NF_model_sex_tissue), sleuth_model_rss(so_TH_model_sex_tissue), sleuth_model_rss(so_NF_t2g_model_sex_tissue), sleuth_model_rss(so_TH_t2g_model_sex_tissue),
                                   sleuth_model_rss(so_NF_model_interaction), sleuth_model_rss(so_TH_model_interaction), sleuth_model_rss(so_NF_t2g_model_interaction), sleuth_model_rss(so_TH_t2g_model_interaction)))

# Post-hoc Target Mapping ----
so_NF_model_sex$target_mapping <- t2g
so_NF_model_sex_tissue$target_mapping <- t2g
so_NF_model_interaction$target_mapping <- t2g
so_TH_model_sex$target_mapping <- t2g
so_TH_model_sex_tissue$target_mapping <- t2g
so_TH_model_interaction$target_mapping <- t2g

# Pseudo-Bulk Analysis ----
## Investigating Histamine and PAR Itch Associated Genes in Sensory Neurons ----
filtered_scaled_transcript_counts(so_NF_t2g_model_interaction, c("Hrh1", "F2rl1"))
# [1] target_id             ens_gene              ext_gene              sample                scaled_reads_per_base tpm                   est_counts
# <0 rows> (or 0-length row.names)

filtered_scaled_transcript_counts(so_NF_model_interaction, c("Hrh1", "F2rl1"))
# [1] target_id             ens_gene              ext_gene              sample                scaled_reads_per_base tpm                   est_counts
# <0 rows> (or 0-length row.names)

# Non-filtered gene expression in sleuth_objects
so_NF_t2g_model_interaction$obs_norm[so_NF_t2g_model_interaction$obs_norm$target_id %in% c('Hrh1', 'F2rl1'),]
so_NF_model_interaction$obs_norm[so_NF_model_interaction$obs_norm$target_id %in% t2g$target_id[t2g$ext_gene %in% c('Hrh1', 'F2rl1')],]

## PLCβ3, TRPV1, and TRPV3 Expression ----
so_NF_t2g_model_interaction$obs_norm[so_NF_t2g_model_interaction$obs_norm$target_id %in% c('Trpv1'),]
so_NF_t2g_model_interaction$obs_norm[so_NF_t2g_model_interaction$obs_norm$target_id %in% c('Trpv3'),]
so_NF_model_interaction$obs_norm[so_NF_model_interaction$obs_norm$target_id %in% t2g$target_id[t2g$ext_gene %in% c('Trpv1')],]
so_NF_model_interaction$obs_norm[so_NF_model_interaction$obs_norm$target_id %in% t2g$target_id[t2g$ext_gene %in% c('Trpv3')],]

# Non-filtered gene expression in sleuth_objects
trpv1_t2g_expression <- so_NF_t2g_model_interaction$obs_norm[so_NF_t2g_model_interaction$obs_norm$target_id %in% c('Trpv1'),]
nrow(trpv1_t2g_expression[trpv1_t2g_expression$scaled_reads_per_base > 5,])
# [1] 12

# Better Transcript Filtering ----
reanalysis_data <- analysis_data
reanalysis_data <- reanalysis_data[-c(2,4),]
reanalysis_data$model_name[1] <- 'NF_model_interaction'; reanalysis_data$model_name[2] <- 'NF_t2g_model_interaction'
reanalysis_data$model_data <- '~sex*tissue'
reanalysis_data$model_parameters[1] <- "filter_fun=function(x){design_filter(metadata, ~sex*tissue, x)}"
reanalysis_data$model_parameters[2] <- "target_map = t2g, aggregation_column = 'ext_gene', gene_mode = TRUE, filter_fun=function(x){design_filter(metadata, ~sex*tissue, x)}"

sleuth_interpret(reanalysis_data, num_core = 1)

# Testing for TRPV1 expression ----
sleuth_test_wt(so_NF_model_interaction)
sleuth_test_wt(so_NF_t2g_model_interaction)

so_NF_model_interaction$target_mapping <- t2g

head(filtered_scaled_transcript_counts(so_NF_model_interaction, "Trpv1"))
head(filtered_scaled_transcript_counts(so_NF_t2g_model_interaction, "Trpv1"))

# Heatmap Visualization and Non-Parametric Experimental Factor Combination Testing ----
# For all the transcripts passing the filter and associated with the gene
heatmap_plot(so_NF_model_interaction, genes = c("Plcb3", "Trpv1"),
             grouping_colours = list(sex = c("F" = "deeppink1", "M" = "dodgerblue2"), tissue = c("NF" = "chartreuse1", "NP" = "blue", "PEP" = "darkorange", "TH" = "gray40")), q_max = FALSE,
             clusterRows = TRUE, clusterColumn = FALSE)

# For all the genes
heatmap_plot(so_NF_t2g_model_interaction, genes = c("Plcb3", "Trpv1"),
             grouping_colours = list(sex = c("F" = "deeppink1", "M" = "dodgerblue2"), tissue = c("NF" = "chartreuse1", "NP" = "blue", "PEP" = "darkorange", "TH" = "gray40")), q_max = FALSE,
             clusterRows = TRUE, clusterColumn = FALSE)

# Non-parametric testing
kw_result <- sleuth_kruskal_wallis(so_NF_model_interaction, "Trpv1")

# Plot the results
ggplot(kw_result, aes(x = est_counts, y = factor_group)) +
  geom_density_ridges2()
