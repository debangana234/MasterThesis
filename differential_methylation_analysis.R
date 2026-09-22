source("package_helpers.R")

run_limma_differential_methylation <- function(methylation_matrix,
                                               metadata,
                                               group_col = "stage",
                                               contrast_strings = c(
                                                 Late_vs_EarlyControl = "stagelate_control - stageearly_control",
                                                 Pre_vs_EarlyControl = "stagePreConverter - stageearly_control",
                                                 Post_vs_LateControl = "stagePostConverter - stagelate_control",
                                                 Post_vs_PreConverter = "stagePostConverter - stagePreConverter"
                                               ),
                                               covariates = c("age", "sex", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"),
                                               match_on_age_sex = FALSE,
                                               sample_col = "Sample_Name",
                                               subject_col = "TREND.ID",
                                               output_dir) {
  install_if_missing(
    packages = c("dplyr", "MatchIt"),
    bioc_packages = c("limma")
  )
  library(dplyr)
  library(MatchIt)
  library(limma)

  analysis_metadata <- metadata

  if (match_on_age_sex) {
    match_df <- analysis_metadata
    match_df$match_group <- as.numeric(match_df[[group_col]] == unique(match_df[[group_col]])[2])
    match_it <- matchit(match_group ~ age + sex, data = match_df, method = "nearest", ratio = 1)
    analysis_metadata <- match.data(match_it)[, colnames(metadata)]
  }

  methylation_matrix <- methylation_matrix[, colnames(methylation_matrix) %in% analysis_metadata[[sample_col]]]
  methylation_matrix <- methylation_matrix[, analysis_metadata[[sample_col]]]

  design_terms <- paste(c(paste0("0 + ", group_col), covariates), collapse = " + ")
  design <- model.matrix(as.formula(paste("~", design_terms)), data = analysis_metadata)

  individual_ids <- factor(analysis_metadata[[subject_col]])
  correlation_estimate <- duplicateCorrelation(methylation_matrix, design = design, block = individual_ids)

  fit <- lmFit(
    methylation_matrix,
    design = design,
    block = individual_ids,
    correlation = correlation_estimate$consensus
  )

  contrast_matrix <- makeContrasts(contrasts = contrast_strings, levels = design)
  fit <- contrasts.fit(fit, contrast_matrix)
  fit <- eBayes(fit)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  results <- list()

  for (contrast_name in colnames(contrast_matrix)) {
    contrast_results <- topTable(fit, coef = contrast_name, adjust.method = "BH", number = Inf)
    contrast_results$CpG_ID <- rownames(contrast_results)
    contrast_results <- contrast_results[, c("CpG_ID", setdiff(colnames(contrast_results), "CpG_ID"))]
    results[[contrast_name]] <- contrast_results

    output_file <- file.path(output_dir, paste0("DMP_", contrast_name, ".csv"))
    write.csv(contrast_results, output_file, row.names = FALSE)
  }

  list(
    metadata = analysis_metadata,
    design = design,
    duplicate_correlation = correlation_estimate,
    fit = fit,
    results = results
  )
}

select_top500_combined_rank_cpgs <- function(dmp_results,
                                             cpg_col = "CpG_ID",
                                             pvalue_col = "adj.P.Val",
                                             logfc_col = "logFC",
                                             top_n = 500,
                                             output_file) {
  dmp_results$rank_adj_pvalues <- rank(dmp_results[[pvalue_col]], ties.method = "average")
  dmp_results$rank_logFC <- rank(-abs(dmp_results[[logfc_col]]), ties.method = "average")
  dmp_results$CombinedRank <- dmp_results$rank_adj_pvalues + dmp_results$rank_logFC

  top_cpgs <- dmp_results[order(dmp_results$CombinedRank), ]
  top_cpgs <- top_cpgs[seq_len(min(top_n, nrow(top_cpgs))), ]

  write.csv(top_cpgs, output_file, row.names = FALSE)

  top_cpgs
}
