# Plotting functions for Gumshoe

# Functions ----
pca_plot <- function(sleuth_obj){
  est_count_matrix <- sleuth:::spread_abundance_by(sleuth_obj$obs_norm_filt, "est_counts", sleuth_obj$sample_to_covariates$sample);
  pca <- prcomp(t(est_count_matrix));
  plot_vector_length <- sqrt((pca$rotation[,1]^2)+(pca$rotation[,2]^2));
  # Highlight the ten longest
  keep <- names(tail(sort(plot_vector_length), 10));
  n <- rownames(pca$rotation);
  pca$rotation[!(rownames(pca$rotation) %in% keep),] <- c(0)
  # Add gene name to the loadings viz of available
  if (!is.null(sleuth_obj$target_mapping)) {
    t2g <- sleuth_obj$target_mapping;
    rownames(t2g) <- t2g$target_id;
    rownames(pca$rotation) <- lapply(n, function(x){ifelse(x %in% keep, paste0(x, "-", t2g[x,2]), ".")})
  }
  autoplot(pca, data=t(est_count_matrix), label=TRUE, shape=FALSE, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3, loadings.label.vjust = 1.7)
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
  if (!"ggrepel" %in% (.packages())){
    cat("Function requires the ggrepel package. If the package is not installed, refer to the gumshoe wiki, otherwise, run the following command to load the package:\n")
    cat("\n")
    cat('library("ggrepel")')
  }
  
  else{
    df <- df[, c(label_column_name, y_axis_column_name, x_axis_column_name)]
    colnames(df) <- c("gene_symbol", "Y_axis", "log2FC")
    
    # Set colours and labels
    df$diff_expressed <- "NO"
    df$diff_expressed[df$log2FC > x_cutoff & df$Y_axis < y_cutoff] <- "UP"
    df$diff_expressed[df$log2FC < -x_cutoff & df$Y_axis < y_cutoff] <- "DOWN"
    
    df$label <- NA
    df$label[df$diff_expressed != "NO"] <- df$gene_symbol[df$diff_expressed != "NO"]
    
    if (any("NO" %in% df$diff_expressed)){
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
    }
    
    else{
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
self_plot_pca <- function(obj, pc_x = 1L, pc_y = 2L, units = 'est_counts',
                          color_by = NULL, graph_name ='Sample Plot Name', scaled = TRUE){
  
  mat <- spread_abundance_by(obj$obs_norm_filt, units,
                             obj$sample_to_covariates$sample)
  
  pca_res <- prcomp(t(mat))
  df <- as.data.frame(pca_res$x[, c(pc_x, pc_y)])
  df$sample <- rownames(df)
  rownames(df) <- NULL
  
  ppv <- plot_pc_variance(obj)
  PCpc <- ppv$data$var
  pc_1 <- paste0('PC1 ', '(', round(PCpc[1], 2), '%)')
  pc_2 <- paste0('PC2 ', '(', round(PCpc[2], 2), '%)')
  
  df <- dplyr::left_join(df, obj$sample_to_covariates,
                         by = 'sample')
  
  ggplot(df, aes(PC1, PC2, colour =  color_by, label = sample)) + 
    theme_bw() +
    theme(legend.position = "none") +
    geom_point(alpha = 0.5) +
    geom_text_repel(show.legend = FALSE) +
    scale_color_manual(name = color_by, values = c("black", "blue", "red", "green")) +
    labs(title = graph_name) +
    xlab(pc_1) +
    ylab(pc_2)
  
  if (scaled){
    ggsave("figure1.png", dpi=300, dev='png', height=PCpc[2] / 5, 
           width=PCpc[1] / 5, 
           units="in")
  }
  
  else{
    ggsave("figure1.png", dpi=300, dev='png', height=4.5, 
           width=6, 
           units="in")
  }
}