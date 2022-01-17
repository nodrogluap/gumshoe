# library(dplyr);  # for "tidyverse" tibble manipulation pipes %>%
# library(ggfortify); # for biplot PCA+loading

std_mean <- function(x){sd(x)/sqrt(length(x))}

sleuth_prep_transcript_dosage <- function(sleuth_obj, metadata, transcript_id){
	# Using the raw observation (scaled) so that the order doesn't matter if we add mutliple dosages in sequence
	transcript_expression_level <- sleuth_obj$obs_raw %>% filter(target_id == transcript_id) %>% mutate(est_counts = est_counts*sleuth_obj$est_counts_sf);
	new_metadata <- metadata %>% left_join(transcript_expression_level[,c("sample","est_counts")], by="sample") %>% rename({{transcript_id}} := est_counts);
	new_model <- as.formula(paste0(deparse(sleuth_obj$full_formula)," + ",transcript_id));
	sleuth_prep(new_metadata, new_model);
}

targetid_u_test_quantitative <- function(sleuth_obj, tid, factor_name, value_threshold){
	d <- sleuth_obj$obs_norm_filt %>% filter(target_id == tid) %>% left_join(sleuth_obj$sample_to_covariates, by="sample");
        # Make a new attribute per sample which is the (hopefully unique) delimited concatenation of the requested factors
        d$kw_group_name <- d[,factor_name] > value_threshold;
        u <- wilcox.test(est_counts ~ kw_group_name, data=d); 
        c(tid,u$p.value,NA);
}

targetid_single_factor_kruskal_wallace_test <- function(sleuth_obj, tid, factor_name){
	d <- sleuth_obj$obs_norm_filt %>% filter(target_id == tid) %>% left_join(sleuth_obj$sample_to_covariates, by="sample"); 
	# Make a new attribute per sample which is the (hopefully unique) delimited concatenation of the requested factors
	d$kw_group_name <- d[,factor_name];
	k <- kruskal.test(est_counts ~ kw_group_name, data=d); 
	c(tid,k$p.value,NA);
}

sleuth_u <- function(sleuth_obj,factor_name,value_threshold, target_ids=unique(sleuth_obj$obs_norm_filt$target_id)){
	U_results_list <- lapply(target_ids, function(target_id){targetid_u_test_quantitative(sleuth_obj, target_id, factor_name, value_threshold)}); 
	U_results_dataframe <- as.data.frame(do.call(rbind, U_results_list));
	U_results_dataframe[,3] <- p.adjust(U_results_dataframe[,2], method="fdr");
	colnames(U_results_dataframe) <- c("target_id", "pval", "qval"); 
	U_results_dataframe$pval <- as.numeric(as.character(U_results_dataframe$pval));
        U_results_dataframe[order(U_results_dataframe$pval),]
}

# If there are mapped read in all replicates, that's considered reliable
sleuth_reliable_target_ids <- function(sleuth_obj, min_count_per_sample = 1){
	sleuth_obj$obs_norm_filt %>% 
  	group_by(target_id) %>% 
  	summarise(min_est_counts = min(est_counts)) %>%
	filter(min_est_counts >= min_count_per_sample) %>%
	select(target_id) %>% pull(1)
}

sleuth_model_rss<- function(sleuth_obj){sum(sleuth_obj$fits$full$summary[,2])}

# This is an ANOVA test without the assumption of equal variance between samples
targetid_one_way_anova_test <- function(sleuth_obj,tid,factor_name){
	d <- sleuth_obj$obs_norm_filt %>% filter(target_id == tid) %>% left_join(sleuth_obj$sample_to_covariates, by="sample"); 
	# Make a new attribute per sample which is the (hopefully unique) delimited concatenation of the requested factors
	d$anova_group_name <- d[,factor_name];
	k <- oneway.test(est_counts ~ anova_group_name, data=d); 
	c(tid,k$p.value,NA);
}

sleuth_oneway_anova <- function(sleuth_obj,factor_name,target_ids=unique(sleuth_obj$obs_norm_filt$target_id)){
	oneway_results_list <- lapply(target_ids, function(target_id){targetid_single_factor_kruskal_wallace_test(sleuth_obj, target_id, factor_name)}); 
	oneway_results_dataframe <- as.data.frame(do.call(rbind, oneway_results_list));
	oneway_results_dataframe[,3] <- p.adjust(oneway_results_dataframe[,2], method="fdr");
	colnames(oneway_results_dataframe) <- c("target_id", "pval", "qval"); 
	oneway_results_dataframe$pval <- as.numeric(as.character(oneway_results_dataframe$pval));
        oneway_results_dataframe[order(oneway_results_dataframe$pval),]
}

sleuth_kruskal_wallace <- function(sleuth_obj,factor_names,target_ids=unique(sleuth_obj$obs_norm_filt$target_id)){
	kw_results_list <- lapply(target_ids, function(target_id){targetid_single_factor_kruskal_wallace_test(sleuth_obj, target_id, factor_names)}); 
	kw_results_dataframe <- as.data.frame(do.call(rbind, kw_results_list));
	kw_results_dataframe[,3] <- p.adjust(kw_results_dataframe[,2], method="fdr");
	colnames(kw_results_dataframe) <- c("target_id", "pval", "qval"); 
	kw_results_dataframe$pval <- as.numeric(as.character(kw_results_dataframe$pval));
        kw_results_dataframe[order(kw_results_dataframe$pval),]
}

targetid_paired_t_test <- function(sleuth_obj, tid, pair_name, factor_name){
	row_by_factor_level <- sleuth_obj$obs_norm_filt %>% filter(target_id == tid) %>% left_join(sleuth_obj$sample_to_covariates, by="sample") %>% group_split(.data[[factor_name]]); 
	# Convoluted way to sort by an unordered factor in tibbles, so we can ensure the pairs are kept together (dplyr's arrange() has conniptions)
	row_by_factor_level[[1]] <- row_by_factor_level[[1]][order(as.character(pull(row_by_factor_level[[1]],pair_name))),];
	row_by_factor_level[[2]] <- row_by_factor_level[[2]][order(as.character(pull(row_by_factor_level[[2]],pair_name))),];
	#print(row_by_factor_level[[1]]);
	#print(row_by_factor_level[[2]]);
	# This t-test assumes log-norml distribution
	t <- t.test(x=log1p(row_by_factor_level[[1]]$est_counts), y=log1p(row_by_factor_level[[2]]$est_counts), paired=TRUE, alternative = "two.sided"); 
	logFCs <- log(row_by_factor_level[[2]]$est_counts/row_by_factor_level[[1]]$est_counts); 
	shapiro <- ifelse(sum(!is.nan(logFCs)) < 3, NA, ifelse(length(unique(logFCs)) < 3, NA, shapiro.test(logFCs)$p.value));
	log_counts <- log1p(row_by_factor_level[[1]]$est_counts);
	# Same order as normal Wald test sleuth_results() output
	c(tid,t$p.value,NA,mean(logFCs),sd(logFCs),mean(log_counts),std_mean(log_counts),shapiro);
	#c(tid,t$p.value,NA,mean(logFCs),sd(logFCs),mean(log_counts),std_mean(log_counts));
}

sleuth_paired_t <- function(sleuth_obj,pair_name,factor_name,target_ids=unique(sleuth_obj$obs_norm_filt$target_id)){
	t_test_results_list <- lapply(target_ids, function(target_id){targetid_paired_t_test(sleuth_obj, target_id, pair_name, factor_name)}); 
	t_test_results_dataframe <- as.data.frame(do.call(rbind, t_test_results_list));
	t_test_results_dataframe[,3] <- p.adjust(t_test_results_dataframe[,2], method="fdr");
	colnames(t_test_results_dataframe) <- c("target_id", "pval", "qval", "b", "se_b", "mean_obs", "var_obs", "Shapiro-Wilks");
	# The following makes it syntactically easier to pull out rows by the callee
	rownames(t_test_results_dataframe) <- t_test_results_dataframe$target_id;
	# Otherwise NAs cause p-value to be considered as a factor column
	t_test_results_dataframe$pval <- as.numeric(as.character(t_test_results_dataframe$pval));
	t_test_results_dataframe[order(t_test_results_dataframe$pval),]
}

fdr_ids <- function(sleuth_obj, test_name, fdr=0.05, test_type="wald"){
	res <- sleuth_results(sleuth_obj, test_name, test_type=test_type);
	res$target_id[which(res$qval <= fdr)];
}

replicate_dotplots <- function(sleuth_obj, design_matrix_combo){
	samples <- sleuth_obj$design_matrix;
	for (column_index in 1:length(design_matrix_combo)){
		# +1 on the index because of the Intercept term in the design matrix
		samples <- samples[which(samples[,column_index+1] == design_matrix_combo[column_index]),];
	} 
	num_replicates <- dim(samples)[1];
	data <- sleuth_obj$obs_norm_filt[which(sleuth_obj$obs_norm_filt$sample %in% rownames(samples)),];
	data$log_est_counts <- log1p(data$est_counts);
	layout(mat=matrix(1:(num_replicates^2),nrow=num_replicates,ncol=num_replicates,byrow=T));
	for (sample1 in rownames(samples)){
		for(sample2 in rownames(samples)){
			plot(data$log_est_counts[which(data$sample == sample1)], data$log_est_counts[which(data$sample == sample2)], xlab=sample1, ylab=sample2, xlim=c(0,12), ylim=c(0,12))
		}
	}
}

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
