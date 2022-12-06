# Plotting functions for Gumshoe

# Functions ----
by_plot <- function(sleuth_obj) {
  est_count_matrix <- sleuth:::spread_abundance_by(sleuth_obj$obs_norm_filt, "est_counts", sleuth_obj$sample_to_covariates$sample)
  pca <- prcomp(t(est_count_matrix))
  plot_vector_length <- sqrt((pca$rotation[, 1]^2) + (pca$rotation[, 2]^2))
  # Highlight the ten longest
  keep <- names(tail(sort(plot_vector_length), 10))
  n <- rownames(pca$rotation)
  pca$rotation[!(rownames(pca$rotation) %in% keep), ] <- c(0)
  # Add gene name to the loadings viz of available
  if (!is.null(sleuth_obj$target_mapping)) {
    t2g <- sleuth_obj$target_mapping
    rownames(t2g) <- t2g$target_id
    rownames(pca$rotation) <- lapply(n, function(x) {
      ifelse(x %in% keep, paste0(x, "-", t2g[x, 2]), ".")
    })
  }
  autoplot(pca, data = t(est_count_matrix), label = TRUE, shape = FALSE, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3, loadings.label.vjust = 1.7)
}

#' Generate a volcano plot from a data frame containing a gene name, q-value, and fold change with user-selected cutoffs and title.
#'
#' @param df A data frame or list containing ensemble transcript id's
#' @param label_column_name The column name of the data point labels. The default is external_gene_name, which is the column name from the output of the ensembl_to_id function.
#' @param x_axis_column_name The column name for the intended x-axis values. The default is "b", which is the name of the column corresponding to the FC after Sleuth analysis.
#' @param y_axis_column_name The column name for the intended y-axis values. The default is "qval", which is the name of the column corresponding to the Benjamini-Hochberg multiple testing corrected p-value after Sleuth analysis.
#' @param x_cutoff X-axis cutoff, which corresponds to the fold change. Default value is 2.
#' @param y_cutoff Y-axis cutoff, which corresponds to the q-value. Default value is -log10(.05).
#' @param graph_name Name of the graph. Default is "Sample Graph Name".
#'
#' @export
#' @examples
#' # Create a volcano plot from a data frame with gene symbols, q-value, and FC.
#' volc_plot(wald_sig_results, graph_name = "Treated v. WT - Sex Factor")
volc_plot <- function(df, label_column_name = "external_gene_name", x_axis_column_name = "b", y_axis_column_name = "qval", x_cutoff = 2, y_cutoff = -log10(.05), graph_name = "Sample Graph Name") {
  if (!"ggrepel" %in% (.packages())) {
    cat("Function requires the ggrepel package. If the package is not installed, refer to the gumshoe wiki, otherwise, run the following command to load the package:\n")
    cat("\n")
    cat('library("ggrepel")')
  } else {
    df <- df[, c(label_column_name, y_axis_column_name, x_axis_column_name)]
    colnames(df) <- c("gene_symbol", "Y_axis", "log2FC")

    # Set colours and labels
    df$diff_expressed <- "NO"
    df$diff_expressed[df$log2FC > x_cutoff & df$Y_axis < y_cutoff] <- "UP"
    df$diff_expressed[df$log2FC < -x_cutoff & df$Y_axis < y_cutoff] <- "DOWN"

    df$label <- NA
    df$label[df$diff_expressed != "NO"] <- df$gene_symbol[df$diff_expressed != "NO"]

    if (any("NO" %in% df$diff_expressed)) {
      ggplot(df, aes(x = log2FC, y = -log10(Y_axis), col = diff_expressed, label = label)) +
        geom_point(alpha = 0.5) +
        theme_minimal() +
        geom_text_repel(show.legend = FALSE, min.segment.length = 0.5) +
        scale_color_manual(name = "Significant", values = c("blue", "black", "red")) +
        geom_vline(xintercept = c(-x_cutoff, x_cutoff), col = "black", linetype = "dashed") +
        geom_hline(yintercept = y_cutoff, col = "black", linetype = "dashed") +
        labs(title = graph_name) +
        xlab(expression(Log[2] ~ FC)) +
        ylab(expression(-Log[10] ~ (Adj. ~ p - Value)))
    } else {
      ggplot(df, aes(x = log2FC, y = -log10(Y_axis), col = diff_expressed, label = label)) +
        geom_point(alpha = 0.5) +
        theme_minimal() +
        geom_text_repel(show.legend = FALSE, min.segment.length = 0.5) +
        scale_color_manual(name = "Significant", values = c("blue", "red")) +
        geom_vline(xintercept = c(-x_cutoff, x_cutoff), col = "black", linetype = "dashed") +
        geom_hline(yintercept = y_cutoff, col = "black", linetype = "dashed") +
        labs(title = graph_name) +
        xlab(expression(Log[2] ~ FC)) +
        ylab(expression(-Log[10] ~ (Adj. ~ p - Value)))
    }
  }
}

#' Spread abundance by a column. Take a data.frame from a sleuth object
#' (e.g. \code{obs_raw}) and cast it into a matrix where the rows are the
#' target_ids and the columns are the sample ids. The values are the variable
#' you are "spreading" on.
#'
#' @param abund the abundance \code{data.frame} from a \code{sleuth} object
#' @param var a character array of length one. The variable for which to get "spread" on (e.g. "est_counts").
#'
#' @export
spread_abundance_by <- function(abund, var, which_order) {
  # var <- lazyeval::lazy(var)
  abund <- data.table::as.data.table(abund)
  var_spread <- data.table::dcast(abund, target_id ~ sample, value.var = var)
  # there is a discrepancy between data table's sorting of character vectors
  # and how tidyr previously (or the order function) sorts character vectors
  # so next step is needed to make sure the order is correct
  var_spread <- var_spread[order(var_spread$target_id), ]
  var_spread <- as.data.frame(var_spread, stringsAsFactors = FALSE)
  rownames(var_spread) <- var_spread$target_id
  var_spread["target_id"] <- NULL
  result <- as.matrix(var_spread)

  result[, which_order, drop = FALSE]
}

#' Generate a PCA plot from a sleuth object with user-selected colour_by and title.
#'
#' @param obj A sleuth object.
#' @param pc_x Integer denoting the principle component to use for the x-axis.
#' @param pc_y Integer denoting the principle component to use for the y-axis.
#' @param units Either 'est_counts' ('scaled_reads_per_base' for gene_mode) or 'tpm'.
#' @param color_by Column name to colour the samples by. Default is NULL.
#' @param graph_name Name of the graph. Default is "Sample Plot Name".
#' @param scaled If the plot should be scaled to the amount of variation. Scale factor is 1:5.
#'
#' @export
#' @examples
#' # Create a PCA plot from a sleuth object.
#' self_plot_pca(so, colour_by = "tissue", graph_name = "All Samples")
self_plot_pca <- function(obj, pc_x = 1L, pc_y = 2L, units = "est_counts",
                          color_by = NULL, graph_name = "Sample Plot Name", scaled = TRUE) {
  mat <- spread_abundance_by(
    obj$obs_norm_filt, units,
    obj$sample_to_covariates$sample
  )

  pca_res <- prcomp(t(mat))
  df <- as.data.frame(pca_res$x[, c(pc_x, pc_y)])
  df$sample <- rownames(df)
  rownames(df) <- NULL

  ppv <- plot_pc_variance(obj)
  PCpc <- ppv$data$var
  pc_1 <- paste0("PC1 ", "(", round(PCpc[1], 2), "%)")
  pc_2 <- paste0("PC2 ", "(", round(PCpc[2], 2), "%)")

  df <- dplyr::left_join(df, obj$sample_to_covariates,
    by = "sample"
  )

  ggplot(df, aes(PC1, PC2, colour = color_by, label = sample)) +
    theme_bw() +
    theme(legend.position = "none") +
    geom_point(alpha = 0.5) +
    geom_text_repel(show.legend = FALSE) +
    scale_color_manual(name = color_by, values = c("black", "blue", "red", "green")) +
    labs(title = graph_name) +
    xlab(pc_1) +
    ylab(pc_2)

  if (scaled) {
    ggsave("figure1.png",
      dpi = 300, dev = "png", height = PCpc[2] / 5,
      width = PCpc[1] / 5,
      units = "in"
    )
  } else {
    ggsave("figure1.png",
      dpi = 300, dev = "png", height = 4.5,
      width = 6,
      units = "in"
    )
  }
}

#' Generate a box plot from kruskal-wallis results with user-selected title.
#'
#' @param kruskal_result A sleuth object.
#' @param graph_name Name of the graph. Default is "Sample Plot Name".
#'
#' @export
#' @examples
#' # Create a box plot from a kruskal-wallis result retrived using sleuth_kruskal_wallis().
#' isoform_boxplot(kw_result, graph_name = "All Isoforms")
isoform_boxplot <- function(kruskal_result, graph_name = "Sample Plot Name") {
  ggplot(kruskal_result, aes(x = target_id, y = est_counts, col = target_id)) +
    theme_minimal() +
    geom_boxplot() +
    guides(color = "none") +
    geom_jitter() +
    scale_fill_brewer(palette = "Dark2") +
    ylab("Estimated Counts") +
    xlab("Target ID") +
    labs(title = graph_name) +
    rotate_x_text(90)
}

#' Generate a isoform level heat map. Can only show a maximum of three genes at the same time.
#'
#' @param sleuth_obj An existing Sleuth object as generated by sleuth_prep() and fit by sleuth_fit(). The transcript name column must be target_id and the gene name column must be 'ext_gene'. Also, the sleuth_obj must have been prepared in non-gene aggregation mode, but with target mapping enabled.
#' @param genes A vector containing gene names whose transcripts should be plotted.
#' @param q_max The maximum q-value that should be used as a cutoff to determine significant transcripts. If set to FALSE, an FDR cutoff will not be applied. Defaults to .05.
#' @param iqf A numeric value between 0 and 1 (inclusive) to denote the quantile of data to display. Defaults to 0, which displays all data.
#' @param test The sleuth test to be analyzed using sleuth_results, can be either "wt" or "lrt". Defaults to "wt".
#' @param sig_marker_values A vector of values that should be used as a cutoff for significance labelling. Defaults to .05, .01, and .001.
#' @param sig_markers Markers of signifigance. Defaults to \*, \*\*, \*\*\*.
#' @param grouping_colours A list of groups and the associated level-colour pairs. See example for syntax. Defaults to empty list.
#' @param boxplot A Boolean value that indicates whether a box plot of scaled transcripts should be displayed.
#'
#' @export
#' @examples
#' # Create a isoform heat map from a sleuth_object.
#' heatmap_plot(so, genes = c("Tp53", "Tnf", "Egfr", "Egr1"), q_max = .01, test = "lrt", grouping_colours = list(sex = c("F" = "thistle1", "M" = "turquoise1"), treatment = c("drug_A" = "indianred1", "drug_B" = "olivedrab1")), boxplot = FALSE)
heatmap_plot <- function(sleuth_obj, genes, q_max = .05, test = "wt", iqf = 0, sig_marker_values = c(.05, .01, .001), sig_markers = c("*", "**", "***"),
                         grouping_colours = list(), plot_title = "Sleuth Heatmap", boxplot = TRUE) {
  # Retrieve the scaled transcript counts for the set of genes passed to the function
  scaled_trancript_counts <- filtered_scaled_transcript_counts(sleuth_obj = sleuth_obj, genes = genes)

  # Check to see if any genes have been retrieved
  if (dim(scaled_trancript_counts)[1] == 0) {
    return("Confirm that the gene name is correct.")
  }

  scaled_trancript_counts$target_id <- paste(scaled_trancript_counts$ext_gene, scaled_trancript_counts$target_id, sep = " - ")

  # Dcast the scaled_trancript_counts to make the columns the sample names, rows the transcript ids, and the cells the est_count values
  scaled_trancript_counts <- dcast(scaled_trancript_counts, target_id ~ sample, value.var = "est_counts")
  rownames(scaled_trancript_counts) <- scaled_trancript_counts$target_id

  # Remove an extra column created from the dcast that contains the transcript ids,
  # take the log2 of the est_counts to plot them in more a representative way,
  # remove the infinite values that occur when taking the log of a zero value,
  # and convert from a df to a matrix.
  scaled_trancript_counts <- scaled_trancript_counts[-1]
  scaled_trancript_counts <- log2(scaled_trancript_counts)
  scaled_trancript_counts[sapply(scaled_trancript_counts, is.infinite)] <- NA
  assign("stc1", scaled_trancript_counts, envir = .GlobalEnv)

  # Calculate the quantile value based on all the est_counts for every transcript for the selected gene with a probability of the user selected iqf parameter
  transcript_mean <- apply(scaled_trancript_counts, 1, function(v) mean(as.numeric(v), na.rm = TRUE))
  quant_val <- quantile(scaled_trancript_counts, probs = iqf, na.rm = TRUE)

  # Remove all the transcripts that do not have an average est_count that is greater than the quant_val
  scaled_trancript_counts <- scaled_trancript_counts[rownames(scaled_trancript_counts) %in% c(names(transcript_mean[transcript_mean > quant_val])), ]
  assign("stc2", scaled_trancript_counts, envir = .GlobalEnv)

  scaled_trancript_count_matrix <- as.matrix(scaled_trancript_counts)

  # Gather the differentially expressed transcripts after running a wt or lrt
  all_results <- sleuth_object_result(sleuth_obj = sleuth_obj, all_data = FALSE, sig_data = FALSE, single_df = TRUE, retrived_from_model = TRUE, test = test)
  all_results <- all_results[all_results$ext_gene %in% genes, ]
  all_results$target_id <- paste(all_results$ext_gene, all_results$target_id, sep = " - ")
  if (is.numeric(q_max)) {
    all_results <- fdr_cutoff(all_results, q_cutoff = q_max)
  }

  # Check to see if any genes have been retrieved
  if (dim(all_results)[1] == 0) {
    return("No transcript of the genes chosen have a significant difference.")
  }

  # Obtain the q-value per transcript.
  transcript_p_val <- dcast(all_results, target_id ~ models, value.var = "qval", fun.aggregate = sum)
  rownames(transcript_p_val) <- transcript_p_val$target_id
  transcript_p_val <- transcript_p_val[-1]

  # Merge the transcript q-values and all transcripts
  transcript_p_val <- merge(transcript_p_val, scaled_trancript_counts, by = "row.names", all.y = TRUE)
  rownames(transcript_p_val) <- transcript_p_val$Row.names
  transcript_p_val <- transcript_p_val[2:(length(unique(all_results$models)) + 1)]
  transcript_p_val[transcript_p_val == 0] <- NA

  # Heatmap Annotation Construction
  ht_opt(legend_border = "black", heatmap_border = TRUE, annotation_border = TRUE)

  if (length(grouping_colours) == 0) {
    ha_top <- HeatmapAnnotation(df = sleuth_obj$sample_to_covariates[, 2:ncol(sleuth_obj$sample_to_covariates)])
  } else {
    ha_top <- HeatmapAnnotation(df = sleuth_obj$sample_to_covariates[, 2:ncol(sleuth_obj$sample_to_covariates)], col = grouping_colours)
  }

  # Generate gene significance heatmap annotation with asterisk.
  qvalue_col_fun <- colorRamp2(c(0, 1), c("white", "white"))
  model_names <- letters[1:length(colnames(transcript_p_val))]
  transcript_p_val[transcript_p_val > max(attributes(qvalue_col_fun)$breaks)] <- NA

  # Determine what symbols should be assigned to the significant markers
  pch <- transcript_p_val
  annotation_symbol <- letters[1:length(sig_marker_values)]
  sig_marker_values <- sort(sig_marker_values)
  for (cutoff in sig_marker_values) {
    pch[pch < as.numeric(cutoff)] <- annotation_symbol[match(cutoff, sig_marker_values)]
  }

  annotations <- na.exclude(as.character(unique(unlist(pch))))
  annotations <- sort(annotations)
  sig_markers <- sort(sig_markers, decreasing = TRUE)
  for (anno in annotations) {
    pch[pch == anno] <- sig_markers[match(anno, annotations)]
  }

  for (model in colnames(transcript_p_val)) {
    if (model != "(Intercept)") {
      letter_index <- match(model, colnames(transcript_p_val))
      letter <- model_names[letter_index]

      ha_temp <- HeatmapAnnotation(model = anno_simple(transcript_p_val[model], col = qvalue_col_fun, pch = as.matrix(pch[model]), pt_size = unit(1, "snpc") * 0.7, width = max_text_width(sig_markers) * 1.2, na_col = "white", border = TRUE), annotation_label = model, which = "row")

      names(ha_temp@anno_list) <- letter
      ha_temp@anno_list[[letter]]@name <- letter

      if (!exists("ha_right")) {
        ha_right <- ha_temp
      } else {
        ha_right <- c(ha_right, ha_temp)
      }
    }
  }

  # Create a legend for the significant markers
  legend_label <- paste("<", sig_marker_values, sep = " ")
  lgd_sig <- Legend(labels = legend_label, border = TRUE, type = "points", pch = sig_markers, title = "q-Value Markers", grid_width = max_text_width(sig_markers) * 1.2, background = qvalue_col_fun(c(sig_marker_values)))

  if (boxplot) {
    boxplot_vals <- apply(scaled_trancript_count_matrix, 1, function(v) median(as.numeric(v), na.rm = TRUE))
    boxplot_col <- colorRamp2(c(min(boxplot_vals), max(boxplot_vals)), c("blue", "red"))
    ha_left <- rowAnnotation("Scaled Transcript Est Counts" = anno_boxplot(scaled_trancript_count_matrix, gp = gpar(fill = boxplot_col(boxplot_vals))))

    ha_c <- Heatmap(scaled_trancript_count_matrix,
      name = "Scaled Counts", rect_gp = gpar(col = "white", lwd = 1), border_gp = gpar(col = "black"), clustering_distance_rows = "pearson", clustering_distance_columns = "pearson",
      row_title = "Transcripts", row_title_rot = 0, row_names_max_width = max_text_width(rownames(scaled_trancript_count_matrix), gp = gpar(fontsize = 12)),
      show_column_names = FALSE,
      na_col = "black",
      top_annotation = ha_top, left_annotation = ha_left, right_annotation = ha_right
    )
  } else {
    ha_c <- Heatmap(scaled_trancript_count_matrix,
      name = "Scaled Counts", rect_gp = gpar(col = "white", lwd = 1), border_gp = gpar(col = "black"), clustering_distance_rows = "pearson", clustering_distance_columns = "pearson",
      row_title = "Transcripts", row_title_rot = 0, row_names_max_width = max_text_width(rownames(scaled_trancript_count_matrix), gp = gpar(fontsize = 12)),
      show_column_names = FALSE,
      na_col = "black",
      top_annotation = ha_top, right_annotation = ha_right
    )
  }
  draw(ha_c, annotation_legend_list = list(lgd_sig), merge_legend = TRUE, column_title = plot_title)
}
