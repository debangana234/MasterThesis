source("Parkinson/clean_scripts/00_package_helpers.R")

apply_blup_horvath_clocks <- function(beta_values,
                                      metadata,
                                      sample_col = "Sample_Name",
                                      age_col = "age",
                                      output_file) {
  install_if_missing(
    packages = c("dplyr", "tibble"),
    bioc_packages = c("methylclock")
  )
  library(dplyr)
  library(tibble)
  library(methylclock)

  methylation_sites_df <- as.data.frame(beta_values)
  methylation_sites_df <- rownames_to_column(methylation_sites_df, "sites")
  colnames(methylation_sites_df) <- gsub("^X", "", colnames(methylation_sites_df))

  age_estimate <- DNAmAge(
    methylation_sites_df,
    age = metadata[[age_col]],
    cell.count = TRUE
  )

  age_df <- as.data.frame(age_estimate)
  age_df <- age_df[, colnames(age_df) %in% c("id", "Horvath", "BLUP")]
  age_df <- merge(age_df, metadata, by.x = "id", by.y = sample_col)

  write.csv(age_df, output_file, row.names = FALSE)

  age_df
}

find_epigenetic_age_acceleration <- function(age_df,
                                             age_col = "age",
                                             clock_cols = c("Horvath", "BLUP"),
                                             celltype_cols = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"),
                                             subject_col = "TREND.ID",
                                             output_file) {
  install_if_missing(packages = c("dplyr", "lmerTest"))
  library(dplyr)
  library(lmerTest)

  result_df <- age_df

  for (clock in clock_cols) {
    covariates <- c(age_col, celltype_cols)

    model_formula <- as.formula(paste(clock, "~", paste(covariates, collapse = " + ")))
    model <- lm(model_formula, data = result_df)
    result_df[[paste0("EAA_", clock)]] <- residuals(model)

    mixed_formula <- as.formula(paste(clock, "~", age_col, "+ (1|", subject_col, ")"))
    mixed_model <- lmer(mixed_formula, data = result_df)
    result_df[[paste0("EAA_", clock, "_mixed")]] <- residuals(mixed_model)
  }

  write.csv(result_df, output_file, row.names = FALSE)

  result_df
}

PCA_biplot <- function(age_df,
                       eaa_col = "EAA_Horvath",
                       sample_col = "Sample_Name",
                       subject_col = "TREND.ID",
                       age_col = "age",
                       group_col = "stage",
                       celltype_cols = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"),
                       stage_levels = c("early_control", "late_control", "PreConverter", "PostConverter"),
                       output_prefix) {
  install_if_missing(packages = c("dplyr", "ggplot2", "ggforce"))
  library(dplyr)
  library(ggplot2)
  library(ggforce)

  plot_data <- age_df
  plot_data[[group_col]] <- factor(plot_data[[group_col]], levels = stage_levels)

  pca_variables <- c(eaa_col, celltype_cols)
  pca_input <- plot_data[, pca_variables]
  complete_rows <- complete.cases(pca_input)
  pca_input <- pca_input[complete_rows, ]
  plot_data <- plot_data[complete_rows, ]

  pca_result <- prcomp(pca_input, center = TRUE, scale. = TRUE)
  variance_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100

  scores <- as.data.frame(pca_result$x[, 1:2])
  scores[[sample_col]] <- plot_data[[sample_col]]
  scores[[subject_col]] <- plot_data[[subject_col]]
  scores[[age_col]] <- plot_data[[age_col]]
  scores[[group_col]] <- plot_data[[group_col]]
  scores[[eaa_col]] <- plot_data[[eaa_col]]

  median_eaa <- median(plot_data[[eaa_col]], na.rm = TRUE)
  scores$EAA_category <- ifelse(scores[[eaa_col]] >= median_eaa, "High", "Low")
  scores$EAA_category <- factor(scores$EAA_category, levels = c("Low", "High"))

  loadings <- as.data.frame(pca_result$rotation[, 1:2])
  loadings$variable <- rownames(loadings)

  score_range_x <- max(scores$PC1, na.rm = TRUE) - min(scores$PC1, na.rm = TRUE)
  score_range_y <- max(scores$PC2, na.rm = TRUE) - min(scores$PC2, na.rm = TRUE)
  loading_range_x <- max(loadings$PC1, na.rm = TRUE) - min(loadings$PC1, na.rm = TRUE)
  loading_range_y <- max(loadings$PC2, na.rm = TRUE) - min(loadings$PC2, na.rm = TRUE)
  arrow_scale <- min(score_range_x / loading_range_x, score_range_y / loading_range_y) * 0.48

  loadings <- loadings %>%
    mutate(arrow_x = PC1 * arrow_scale,
      arrow_y = PC2 * arrow_scale,
      label_x = arrow_x * 1.12,
      label_y = arrow_y * 1.12)

  variance_table <- data.frame(
    PC = paste0("PC", seq_along(variance_explained)),
    variance_percent = variance_explained
  )

  pca_plot <- ggplot(scores, aes(x = PC1, y = PC2)) +
    geom_mark_ellipse(aes(fill = EAA_category, group = EAA_category), alpha = 0.12, color = NA, show.legend = FALSE) +
    geom_hline(yintercept = 0, color = "grey82", linewidth = 0.45) +
    geom_vline(xintercept = 0, color = "grey82", linewidth = 0.45) +
    geom_point(aes(fill = .data[[group_col]], shape = EAA_category), color = "white", size = 3, stroke = 0.55, alpha = 0.94) +
    geom_segment(data = loadings,
      aes(x = 0, y = 0, xend = arrow_x, yend = arrow_y),
      inherit.aes = FALSE,
      arrow = arrow(length = unit(0.20, "cm")),
      color = "#27343B",
      linewidth = 0.75) +
    geom_text(data = loadings,
      aes(x = label_x, y = label_y, label = variable),
      inherit.aes = FALSE,
      size = 3.3,
      color = "#27343B") +
    theme_minimal() +
    labs(x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
      y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
      fill = "Stage",
      shape = "EAA")

  write.csv(scores, paste0(output_prefix, "_scores.csv"), row.names = FALSE)
  write.csv(loadings, paste0(output_prefix, "_loadings.csv"), row.names = FALSE)
  write.csv(variance_table, paste0(output_prefix, "_variance_explained.csv"), row.names = FALSE)
  ggsave(paste0(output_prefix, ".png"), plot = pca_plot, width = 8, height = 6, dpi = 300)
  ggsave(paste0(output_prefix, ".pdf"), plot = pca_plot, width = 8, height = 6)

  list(pca_result = pca_result,
    scores = scores,
    loadings = loadings,
    variance_explained = variance_table,
    plot = pca_plot)}

plot_epigenetic_age_boxplot <- function(age_df,
                                        eaa_col = "EAA_BLUP",
                                        group_col = "stage",
                                        stage_levels = c("early_control", "late_control", "PreConverter", "PostConverter"),
                                        output_file) {
  install_if_missing(packages = c("ggplot2"))
  library(ggplot2)

  age_df[[group_col]] <- factor(age_df[[group_col]], levels = stage_levels)

  plot <- ggplot(age_df, aes(x = .data[[group_col]], y = .data[[eaa_col]], fill = .data[[group_col]])) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_point(color = "blue", alpha = 0.7) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(x = "", y = eaa_col)

  ggsave(output_file, plot = plot, width = 7, height = 5, dpi = 300)

  plot
}
