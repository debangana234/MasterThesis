source("package_helpers.R")

plot_cell_type_proportions <- function(celltype_df,
                                       metadata,
                                       sample_col = "Sample_Name",
                                       group_col = "stage",
                                       cell_types = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"),
                                       stage_levels = c("early_control", "late_control", "PreConverter", "PostConverter"),
                                       output_file = NULL) {
  install_if_missing(packages = c("dplyr", "tidyr", "ggplot2", "rstatix"))
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rstatix)

  celltypemetadata <- merge(metadata, celltype_df, by = sample_col)
  celltypemetadata[[group_col]] <- factor(celltypemetadata[[group_col]], levels = stage_levels)

  celltypes_long <- celltypemetadata |>
    pivot_longer(
      cols = all_of(cell_types),
      names_to = "CellType",
      values_to = "Proportion"
    )

  plot <- ggplot(celltypes_long, aes(x = .data[[group_col]], y = Proportion, fill = .data[[group_col]])) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(color = "blue", width = 0.08, alpha = 0.7, size = 1.5) +
    facet_wrap(~CellType, scales = "free_y") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    labs(x = "", y = "Cell type proportion")

  pairwise_results <- data.frame()

  for (cell_type in unique(celltypes_long$CellType)) {
    subset_data <- celltypes_long[celltypes_long$CellType == cell_type, ]

    test_results <- pairwise.wilcox.test(
      subset_data$Proportion,
      subset_data[[group_col]],
      p.adjust.method = "BH"
    )

    p_values <- as.data.frame(as.table(test_results$p.value))
    p_values <- p_values[!is.na(p_values$Freq), ]
    p_values$CellType <- cell_type
    colnames(p_values) <- c("Var1", "Var2", "p_value", "CellType")

    pairwise_results <- rbind(pairwise_results, p_values)
  }

  if (!is.null(output_file)) {
    ggsave(output_file, plot = plot, width = 10, height = 7, dpi = 300)
  }

  list(
    merged_data = celltypemetadata,
    long_data = celltypes_long,
    pairwise_results = pairwise_results,
    plot = plot
  )
}

model_longitudinal_proportion_variation <- function(celltype_df,
                                                    metadata,
                                                    sample_col = "Sample_Name",
                                                    subject_col = "TREND.ID",
                                                    time_col = "diagnosis_time",
                                                    group_col = "stage",
                                                    age_col = "age",
                                                    sex_col = "sex",
                                                    cell_types = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"),
                                                    output_file = NULL) {
  install_if_missing(packages = c("dplyr", "lmerTest"))
  library(dplyr)
  library(lmerTest)

  df_merge <- merge(metadata, celltype_df, by = sample_col)
  all_results <- data.frame()

  for (cell in cell_types) {
    model_formula <- as.formula(
      paste(cell, "~", time_col, "+", age_col, "+", sex_col, "+", group_col, "+ (1|", subject_col, ")")
    )

    fit <- lmer(model_formula, data = df_merge)
    fit_summary <- as.data.frame(summary(fit)$coefficients)
    fit_summary$covariate <- rownames(fit_summary)
    fit_summary$cell <- cell
    rownames(fit_summary) <- NULL

    all_results <- rbind(all_results, fit_summary)
  }

  if (!is.null(output_file)) {
    write.csv(all_results, output_file, row.names = FALSE)
  }

  all_results
}

changing_immune_proportions_diagnosis_time <- function(celltype_df,
                                                       metadata,
                                                       sample_col = "Sample_Name",
                                                       subject_col = "TREND.ID",
                                                       age_col = "age",
                                                       diagnosis_age_col = "diagnosis_age",
                                                       group_col = "stage",
                                                       converter_stages = c("PreConverter", "PostConverter"),
                                                       cell_types = c("Neu", "NK", "Mono"),
                                                       cell_type_labels = c(Neu = "Neutrophils", NK = "NK cells", Mono = "Monocytes"),
                                                       output_prefix) {
  install_if_missing(packages = c("dplyr", "tidyr", "ggplot2"))
  library(dplyr)
  library(tidyr)
  library(ggplot2)

  celltypemetadata <- merge(metadata, celltype_df, by = sample_col)

  converter_data <- celltypemetadata %>%
    filter(.data[[group_col]] %in% converter_stages)

  converter_data$Subject <- converter_data[[subject_col]]
  converter_data$Stage <- factor(converter_data[[group_col]], levels = converter_stages)
  converter_data$YearsRelativeDiagnosis <- as.numeric(converter_data[[age_col]]) - as.numeric(converter_data[[diagnosis_age_col]])

  converter_long <- converter_data %>%
    pivot_longer(
      cols = all_of(cell_types),
      names_to = "CellType",
      values_to = "Proportion"
    )

  converter_long$CellTypeLabel <- cell_type_labels[converter_long$CellType]
  converter_long$CellTypeLabel <- factor(converter_long$CellTypeLabel, levels = cell_type_labels[cell_types])

  slope_results <- data.frame()

  for (cell_type in levels(converter_long$CellTypeLabel)) {
    cell_data <- converter_long[converter_long$CellTypeLabel == cell_type, ]
    slope_model <- lm(Proportion ~ YearsRelativeDiagnosis, data = cell_data)
    slope_table <- as.data.frame(summary(slope_model)$coefficients)
    slope_table$term <- rownames(slope_table)
    slope_table$CellType <- cell_type
    rownames(slope_table) <- NULL
    slope_results <- rbind(slope_results, slope_table)
  }

  slope_labels <- slope_results %>%
    filter(term == "YearsRelativeDiagnosis") %>%
    select(CellType, Estimate) %>%
    mutate(
      label = paste0("slope/year = ", ifelse(Estimate >= 0, "+", ""), round(Estimate, 4)),
      YearsRelativeDiagnosis = -Inf,
      Proportion = Inf
    )

  progression_plot <- ggplot(converter_long, aes(x = YearsRelativeDiagnosis, y = Proportion)) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#ead7ff", alpha = 0.35) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.7) +
    geom_line(aes(group = Subject, color = Subject), alpha = 0.6, linewidth = 0.7) +
    geom_point(aes(fill = Stage), shape = 21, size = 2.4, color = "grey35", alpha = 0.95) +
    geom_smooth(method = "lm", se = FALSE, color = "#b00000", linetype = "dashed", linewidth = 1.1) +
    geom_text(
      data = slope_labels,
      aes(x = YearsRelativeDiagnosis, y = Proportion, label = label),
      inherit.aes = FALSE,
      color = "#9b2d27",
      fontface = "bold",
      size = 3.5,
      hjust = -0.1,
      vjust = 1.5
    ) +
    facet_wrap(~CellTypeLabel, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = c("PreConverter" = "#e85a91", "PostConverter" = "#d85c3a")) +
    guides(color = "none") +
    theme_minimal() +
    theme(
      strip.text = element_text(face = "bold", size = 13),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    ) +
    labs(
      title = "Cell Composition Relative to Diagnosis in Converters",
      subtitle = "Converters aligned to diagnosis age.",
      x = "Years relative to diagnosis",
      y = "Cell type proportion",
      fill = "Stage"
    )

  dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
  write.csv(converter_long, paste0(output_prefix, "_diagnosis_aligned_celltypes.csv"), row.names = FALSE)
  write.csv(slope_results, paste0(output_prefix, "_diagnosis_time_slopes.csv"), row.names = FALSE)
  ggsave(paste0(output_prefix, "_diagnosis_aligned_plot.png"), progression_plot, width = 12, height = 6, dpi = 300)
  ggsave(paste0(output_prefix, "_diagnosis_aligned_plot.pdf"), progression_plot, width = 12, height = 6)

  list(
    plot_data = converter_long,
    slope_results = slope_results,
    plot = progression_plot
  )
}
