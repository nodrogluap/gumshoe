# Additional pre-processing functions to complement Gumshoe for Sleuth analysis.

# Required libraries ----
# library(tidyverse)
# library(sleuth)
# library(biomaRt)


# Parameters ----
# Generate the necessary models for sleuth to run. There are five in total:
# model_nominal, model_quant_linear/log/exp, and model_two_binary_factors.

# Define the metadata files, associated models, and the corresponding model
# formula to be used. The `I` before the mathematical expressions makes the
# expression uninterpreted until Sleuth processes the model. NOTE: The order
# of the model data must be align with the model names.
# metadata_names <- c()
# model_names <- c()
# model_data <- c()


# Analysis data frame creation and metadata manipulation ----
# The analysis data file contains the metadata file name, the metadata file, the
# model names, and the model formula.
# analysis_data <- data.frame(metadata_name = metadata_names, metadata_file = tibble(list(metadata)), model_name = model_names, model_data = model_data)
# View(analysis_data)


# Functions ----
#' Automate the running of sleuth_prep and sleuth_fit on all provided metadata files and the associated models
#'
#' @param data A data frame containing metadata file name, metadata tibble file, model names, and model formulae
#' @param num_core An integer for the number of cores to be used for sleuth_prep
#' @param ... Can add in other parameters that will be passed to sleuth_prep
#'
#' @return Sleuth object named using the associated metadata file name and processed with sleuth_prep() and sleuth_fit()
#' @export
#' @examples
#' # Given a sample data frame containing a metadata_file_name, metadata_file, model_name, and model_data (formulae),
#' # create a sleuth object that has been processed using sleuth_prep() followed by sleuth_fit()
#' sleuth_interpret(analysis_data, num_core = 3)
sleuth_interpret <- function(data, num_core = 1, ...) {
  extra_opts <- list(...)

  for (metadata_file_number in 1:length(data$metadata_name)) {
    metadata_file <- data.frame(data[metadata_file_number, 2])
    metadata_model_names <- unlist(strsplit(data$model_name[metadata_file_number], ","))
    metadata_model_formula <- unlist(strsplit(data$model_data[metadata_file_number], ","))

    for (metadata_model_number in 1:length(metadata_model_names)) {
      so_holder_variable <- sleuth_prep(
        sample_to_covariates = metadata_file,
        full_model = as.formula(metadata_model_formula[metadata_model_number]),
        num_cores = num_core,
        ... = extra_opts)
      }

      so_holder_variable <- sleuth_fit(obj = so_holder_variable)

      sleuth_obj_name <- paste("so", metadata_model_names[metadata_model_number], sep = "_")

      assign(sleuth_obj_name, so_holder_variable, envir = .GlobalEnv)
    }
  }

#' Automated function to run all possible Wald tests on a given Sleuth object and the fitted model
#'
#' @param sleuth_obj An existing Sleuth object as generated by sleuth_prep() and fit by sleuth_fit()
#'
#' @return Results of the Wald test assigned to the original Sleuth object input
#' @export
#' @examples
#' # Given a Sleuth object, run a Wald test on all possible models as determined on the formula used to run sleuth_prep()
#' sleuth_test_wt(so)
sleuth_test_wt <- function(sleuth_obj) {
  sleuth_obj_name <- deparse(substitute(sleuth_obj))

  coeff <- colnames(sleuth_obj$design_matrix)

  for (model in coeff) {
    if (model != "(Intercept)") {
      sleuth_obj <- sleuth_wt(sleuth_obj, model)

      assign(sleuth_obj_name, sleuth_obj, envir = .GlobalEnv)
    }
  }
}

#' Automated function to run all possible Likewise Ratio tests on a given Sleuth object with models derived from the original formula used to run sleuth_fit()
#'
#' @param sleuth_obj An existing Sleuth object as generated by sleuth_prep() and fit by sleuth_fit()
#'
#' @return Results of the Likewise Ratio test assigned to the original Sleuth object input
#' @export
#' @examples
#' # Given a Sleuth object, run a Likewise Ratio test on all possible models derived from the original formula used to run sleuth_prep()
#' sleuth_test_lrt(so)
sleuth_test_lrt <- function(sleuth_obj) {
  sleuth_obj_name <- deparse(substitute(sleuth_obj))

  formula_variables <- sleuth_obj$full_formula
  formula_variables <- all.vars(formula_variables)

  for (model in formula_variables) {
    remaining_variables <- paste(formula_variables[-(match(model, formula_variables))], collapse = "_")
    remaining_variables <- paste("no", remaining_variables, sep = "_")

    model <- paste("~", model, sep = "")

    sleuth_obj <- sleuth_fit(sleuth_obj, as.formula(model), fit_name = remaining_variables)
    sleuth_obj <- sleuth_lrt(sleuth_obj, remaining_variables, "full")
    assign(sleuth_obj_name, sleuth_obj, envir = .GlobalEnv)
  }
}

#' Retrieve all target id's for a given data frame with a FDR less than or equal to the selected cutoff. Default cutoff value is .05.
#'
#' @param sleuth_res A data frame generated following running sleuth_results
#' @param q_cutoff The selected FDR cutoff value. Default is .05.
#' @param q_equal Should the returned results by less than or equal to the q_cutoff. Default is false, which then returns all results less than but not equal to the q_cutoff.
#'
#' @return A list of target id's.
#' @export
#' @examples
#' # Given a Sleuth object, return a data frame of all target id's with an FDR less than the cutoff value
#' fdr_cutoff(wald_test_results)
fdr_cutoff <- function(sleuth_res, q_cutoff = .05, q_equal = FALSE) {
  if (q_equal) {
    sleuth_res[which(sleuth_res$qval <= q_cutoff), ]
  } else {
    sleuth_res[which(sleuth_res$qval < q_cutoff), ]
  }
}

#' Retrieve all results from a given Sleuth object following statistical testing
#'
#' @param sleuth_obj An existing Sleuth object as generated by sleuth_prep() and fit by sleuth_fit()
#' @param all_data A Boolean value that indicates whether all results of the Wald test should be assigned to a variable. Automatic assignment to a variable in the global environment. Defaults to FALSE
#' @param sig_data A Boolean value that indicates whether all the significant results of the Wald test coupled with fdr_cutoff should be assigned to a variable. Automatic assignment to a variable in the global environment. Defaults to TRUE.
#' @param single_df A Boolean value that indicates whether all the results of Wald test should be merged into a single data frame and returned. Returned data must be manually assigned to variable when calling with single_df = TRUE. Defaults to FALSE.
#' @param retrived_from_model A Boolean value that indicates whether the data frame returned should include an additional column that includes the model that the data is from. This parameter is unique to single_df as both all_data and sig_data automatically assign data to variables in the global environment. Defaults to FALSE.
#' @param q_max The FDR cutoff value that is passed to the fdr_cutoff function embedded within sleuth_object_result. Defaults to .05.
#' @param test The sleuth test to be analyzed using sleuth_results, can be either "wt" or "lrt". Defaults to "wt".
#'
#' @returns Either automatically generated variables (all_data and sig_data) containing resulting data from the Wald or Likewise Ratio tests with or without an FDR cutoff applied or all test result data in a manually assigned single data frame
#' @export
#' @examples
#' # Given a Sleuth object and Boolean values of tests to run, perform the requested tests on all possible models for the given Sleuth object
#' sleuth_object_result(so, q_max = .01, test = "lrt")
sleuth_object_result <- function(sleuth_obj, all_data = FALSE, sig_data = TRUE, single_df = FALSE, retrived_from_model = FALSE, q_max = .05, test = "wt") {
  coeff <- colnames(sleuth_obj$design_matrix)
  sleuth_obj_name <- deparse(substitute(sleuth_obj))
  sleuth_obj_name <- unlist(strsplit(sleuth_obj_name, "_"))
  sleuth_obj_name <- paste(sleuth_obj_name[c(2:length(sleuth_obj_name))], collapse = "_")

  all_results_single_df <- list()

  if (test == "wt") {
    for (model in coeff) {
      if (model != "(Intercept)") {
        wald_result <- sleuth_results(sleuth_obj, model, test_type = "wt")

        # If the model is an interaction term, remove the colon
        model <- gsub(":", "_", model)

        wald_model_name <- paste("wald", sleuth_obj_name, sep = "_")
        wald_model_name <- paste(wald_model_name, model, sep = "_")

        if (all_data) {
          assign(wald_model_name, wald_result, envir = .GlobalEnv)
        }
        if (sig_data) {
          sig_target_ids <- fdr_cutoff(wald_result, q_cutoff = q_max)
          sig_model_name <- paste("sig", wald_model_name, sep = "_")
          assign(sig_model_name, sig_target_ids, envir = .GlobalEnv)
        }
        if (single_df) {
          if (retrived_from_model) {
            wald_result$models <- model
            all_results_single_df <- rbind(all_results_single_df, wald_result)
          } else {
            all_results_single_df <- rbind(all_results_single_df, wald_result)
          }
        }
      }
    }
  }

  if (test == "lrt") {
    for (model in names(sleuth_obj$tests[["lrt"]])) {
      lrt_result <- sleuth_results(sleuth_obj, model, test_type = "lrt")

      # If the model is an interaction term, remove the colon
      model <- gsub(":", "_", model)

      lrt_model_name <- paste("lrt", sleuth_obj_name, sep = "_")
      lrt_model_name <- paste(lrt_model_name, model, sep = "_")

      if (all_data) {
        assign(lrt_model_name, lrt_result, envir = .GlobalEnv)
      }
      if (sig_data) {
        sig_target_ids <- fdr_cutoff(lrt_result, q_cutoff = q_max)
        sig_model_name <- paste("sig", lrt_model_name, sep = "_")
        assign(sig_model_name, sig_target_ids, envir = .GlobalEnv)
      }
      if (single_df) {
        if (retrived_from_model) {
          lrt_result$models <- model
          all_results_single_df <- rbind(all_results_single_df, lrt_result)
        } else {
          all_results_single_df <- rbind(all_results_single_df, lrt_result)
        }
      }
    }
  }
  if (single_df) {
    all_results_single_df
  }
}


#' Convert a given list of ensemble transcript id's into \code{ensembl_transcript_id}'s, \code{external_gene_name}'s, \code{ensembl_gene_id}'s, and the original transcript.id
#'
#' @param sig_results A data frame or list containing ensemble transcript id's
#' @param entire_gene_name Set to FALSE if the non-abbreviated gene name is required. Default value is TRUE.
#'
#' @return Converts the original data frame or list into a data frame with columns corresponding to \code{external_gene_name}'s, \code{description}'s, and the original \code{ensembl_transcript_id_version}.
#' @export
#' @examples
#' # Convert a data frame of ensemble transcript id's
#' ensemble_to_id(wald_sig_results)
ensembl_to_id <- function(sig_results, entire_gene_name = TRUE) {
  if (!exists("mart")) {
    cat("Function requires a BioMart database and dataset. Refer to the gumshoe wiki on details to install biomaRt.\n")
    cat("\n")
    cat("If biomaRt is installed, run the following command to create the required `mart` object for mice:\n")
    cat("\n")
    cat('mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")\n')
    cat("\n")
    cat("If a mart object has already been assigned, please confirm that the object is named `mart`")
  } else {
    ensembl_convert <- getBM(
      attributes = c("ensembl_transcript_id_version", "external_gene_name", "description"),
      filters = "ensembl_transcript_id_version",
      values = sig_results$target_id,
      mart = mart
    )

    sig_results <- merge(sig_results, ensembl_convert, by.x = "target_id", by.y = "ensembl_transcript_id_version")

    sig_results <- sig_results[, c(length(sig_results) - 1, 1:(length(sig_results) - 2), length(sig_results))]

    if (entire_gene_name) {
      sig_results <- sig_results[]
    } else {
      sig_results <- sig_results[, c(1:12)]
    }
    assign("annotated_sig_results", sig_results, envir = .GlobalEnv)
  }
}
