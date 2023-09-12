# Libraries ----
library(tidyverse)
library(sleuth)
library(biomaRt)
library(gumshoe)
library(ggrepel)
setwd("~/Documents/gumshoe/tutorial")

# Metadata Files ----
## Load Metadata File ----
default_metadata <- read.table("metadata.txt", header = TRUE)

## Generate Metadata File ----
# Create a df that is 4 columns by 24 rows and rename the columns
metadata <- data.frame(matrix(NA, nrow = 24, ncol = 4))
colnames(metadata) <- c("sample", "path", "sex", "tissue")

# Assign sample names. The naming format is sex, tissue type, and replicate number.
metadata[1] <- c("female_NF.0", "female_NF.1", "female_NF.2", "female_NP.0", "female_NP.1", "female_NP.2", "female_PEP.0", "female_PEP.1", "female_PEP.2", "female_TH.0", "female_TH.1", "female_TH.2", "male_NF.0", "male_NF.1", "male_NF.2", "male_NP.0", "male_NP.1", "male_NP.2", "male_PEP.0", "male_PEP.1", "male_PEP.2", "male_TH.0", "male_TH.1", "male_TH.2")

# Assign the file path
# NOTE: This line of code assumes that your data is contained in the working directory in a subfolder named "output".
for (cell in 1:24){
  metadata[cell, 2] <- paste("output/", metadata[cell, 1], sep = "")
  metadata[cell, 2] <- paste(metadata[cell, 2], ".kallisto", sep = "")
}

# Assign sex
metadata[1:12, 3] <- "F"
metadata[13:24, 3] <- "M"

# Assign the tissue type
for (cell in 1:24){
  metadata[cell, 4] <- unlist(strsplit(metadata[cell, 1], "_"))[2]
  metadata[cell, 4] <- gsub("\\.[0-9]", "", metadata[cell, 4])
}

## Compare Contained and Generated Metadata ----
metadata == default_metadata



# Gumshoe Data Frame Creation ----
## Re-level Metadata ----
# Convert the tissue type column to a factor and, based on how R works, NF will be the first factor level.
metadata$tissue <- as.factor(metadata$tissue)

# Create a new variable that will contain a relevled version of the metadata
metadata_releveled <- metadata
metadata_releveled$tissue <- relevel(metadata_releveled$tissue, ref = "TH")


## Analysis Information ----
# Define the metadata file name(s), associated model name(s), and the corresponding model(s) to be used. The `I` before the mathematical expressions informs R that the expression should remain uninterpreted until Sleuth processes the model. 
# NOTE: The order of the model data must be align with the model names.
metadata_names <- c("metadata_nominal_NF_baseline",
                    "metadata_nominal_TH_baseline")
model_names <- c("NF_model_sex,NF_model_sex_tissue,NF_model_interaction",
                 "TH_model_sex,TH_model_sex_tissue,TH_model_interaction")
model_data <- c("~sex, ~sex + tissue, ~sex*tissue",
                "~sex, ~sex + tissue, ~sex*tissue")


## Data Frame Creation ----
# Take the metadata names, model names, and model data information and combine it into a single data frame.
analysis_data <- data.frame(metadata_name = metadata_names, metadata_file = tibble(list(metadata)), model_name = model_names, model_data = model_data)
analysis_data[[2]][[2]] <- metadata_releveled


## Verify the Metadata in the Created Data Frame ----
levels(analysis_data[[2]][[1]]$tissue)
# "NF"  "NP"  "PEP" "TH" 
levels(analysis_data[[2]][[2]]$tissue)
# "TH"  "NF"  "NP"  "PEP"



# Gumshoe Function Calls ----
sleuth_interpret(analysis_data)

# Run a Wald test on each sleuth_object
sleuth_test_wt(so_NF_model_sex)
sleuth_test_wt(so_NF_model_sex_tissue)
sleuth_test_wt(so_NF_model_interaction)
sleuth_test_wt(so_TH_model_sex)
sleuth_test_wt(so_TH_model_sex_tissue)
sleuth_test_wt(so_TH_model_interaction)

# Obtain the results from the sleuth_object
sleuth_object_result(so_NF_model_sex, all_data = FALSE, q_max = 0.01)
sleuth_object_result(so_NF_model_sex_tissue, all_data = FALSE, q_max = 0.01)
sleuth_object_result(so_NF_model_interaction, all_data = FALSE, q_max = 0.01)
sleuth_object_result(so_TH_model_sex, all_data = FALSE, q_max = 0.01)
sleuth_object_result(so_TH_model_sex_tissue, all_data = FALSE, q_max = 0.01)
sleuth_object_result(so_TH_model_interaction, all_data = FALSE, q_max = 0.01)

# Save all the results.
# Save the NF baseline results as .txt in the current directory
write.table(sig_wald_NF_model_interaction_sexM, row.names = FALSE, sep = "\t", quote = FALSE, "NF_interaction_model_sexM.txt")
write.table(sig_wald_NF_model_interaction_sexM_tissueNP, row.names = FALSE, sep = "\t", quote = FALSE, "NF_interaction_model_sexM_NP.txt")
write.table(sig_wald_NF_model_interaction_sexM_tissuePEP, row.names = FALSE, sep = "\t", quote = FALSE, "NF_interaction_model_sexM_PEP.txt")
write.table(sig_wald_NF_model_interaction_sexM_tissueTH, row.names = FALSE, sep = "\t", quote = FALSE, "NF_interaction_model_sexM_TH.txt")
write.table(sig_wald_NF_model_interaction_tissueNP, row.names = FALSE, sep = "\t", quote = FALSE, "NF_interaction_model_NP.txt")
write.table(sig_wald_NF_model_interaction_tissuePEP, row.names = FALSE, sep = "\t", quote = FALSE, "NF_interaction_model_PEP.txt")
write.table(sig_wald_NF_model_interaction_tissueTH, row.names = FALSE, sep = "\t", quote = FALSE, "NF_interaction_model_TH.txt")
write.table(sig_wald_NF_model_sex_sexM, row.names = FALSE, sep = "\t", quote = FALSE, "NF_sex_model_sexM.txt")
write.table(sig_wald_NF_model_sex_tissue_sexM, row.names = FALSE, sep = "\t", quote = FALSE, "NF_sex_tissue_model_sexM.txt")
write.table(sig_wald_NF_model_sex_tissue_tissueNP, row.names = FALSE, sep = "\t", quote = FALSE, "NF_sex_tissue_model_NP.txt")
write.table(sig_wald_NF_model_sex_tissue_tissuePEP, row.names = FALSE, sep = "\t", quote = FALSE, "NF_sex_tissue_model_PEP.txt")
write.table(sig_wald_NF_model_sex_tissue_tissueTH, row.names = FALSE, sep = "\t", quote = FALSE, "NF_sex_tissue_model_TH.txt")

# Save the TH baseline results as .txt in the current directory
write.table(sig_wald_TH_model_interaction_sexM, row.names = FALSE, sep = "\t", quote = FALSE, "TH_interaction_model_sexM.txt")
write.table(sig_wald_TH_model_interaction_sexM_tissueNF, row.names = FALSE, sep = "\t", quote = FALSE, "TH_interaction_model_sexM_NF.txt")
write.table(sig_wald_TH_model_interaction_sexM_tissueNP, row.names = FALSE, sep = "\t", quote = FALSE, "TH_interaction_model_sexM_NP.txt")
write.table(sig_wald_TH_model_interaction_sexM_tissuePEP, row.names = FALSE, sep = "\t", quote = FALSE, "TH_interaction_model_sexM_PEP.txt")
write.table(sig_wald_TH_model_interaction_tissueNF, row.names = FALSE, sep = "\t", quote = FALSE, "TH_interaction_model_NF.txt")
write.table(sig_wald_TH_model_interaction_tissueNP, row.names = FALSE, sep = "\t", quote = FALSE, "TH_interaction_model_NP.txt")
write.table(sig_wald_TH_model_interaction_tissuePEP, row.names = FALSE, sep = "\t", quote = FALSE, "TH_interaction_model_PEP.txt")
write.table(sig_wald_TH_model_sex_sexM, row.names = FALSE, sep = "\t", quote = FALSE, "TH_sex_model_sexM.txt")
write.table(sig_wald_TH_model_sex_tissue_sexM, row.names = FALSE, sep = "\t", quote = FALSE, "TH_sex_tissue_model_sexM.txt")
write.table(sig_wald_TH_model_sex_tissue_tissueNF, row.names = FALSE, sep = "\t", quote = FALSE, "TH_sex_tissue_model_NF.txt")
write.table(sig_wald_TH_model_sex_tissue_tissueNP, row.names = FALSE, sep = "\t", quote = FALSE, "TH_sex_tissue_model_NP.txt")
write.table(sig_wald_TH_model_sex_tissue_tissuePEP, row.names = FALSE, sep = "\t", quote = FALSE, "TH_sex_tissue_model_PEP.txt")

# Create a mart object to retrive the required information
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Call the ensembl_to_id function on a data frame containing results from the NF baseline metadata file with
# the model ~sex
ensembl_to_id(sig_wald_NF_model_sex_sexM)

# Create a volcano plot of the annotated data
volc_plot(annotated_sig_results, graph_name = "Sample Graph - NF Baseline, ~sex, sexM")

