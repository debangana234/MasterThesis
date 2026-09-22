source("package_helpers.R")

plot_volcano_and_dmp_boxplots <- function(dmp_results,
                                          beta_values,
                                          metadata,
                                          sample_col = "Sample_Name",
                                          group_col = "stage",
                                          cpg_col = "CpG_ID",
                                          pvalue_col = "adj.P.Val",
                                          logfc_col = "logFC",
                                          selected_cpgs,
                                          output_dir = NULL) {
  install_if_missing(packages = c("dplyr", "tidyr", "ggplot2", "tibble"))
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)

  dmp_results$minus_log10_p <- -log10(dmp_results[[pvalue_col]])
  dmp_results$significant <- dmp_results[[pvalue_col]] < 0.05

  volcano_plot <- ggplot(dmp_results, aes(x = .data[[logfc_col]], y = minus_log10_p, color = significant)) +
    geom_point(alpha = 0.7, size = 1.5) +
    theme_minimal() +
    labs(x = "logFC", y = "-log10 adjusted P value")

  boxplot_list <- list()

  for (cpg_site in selected_cpgs) {
    methylation_significant_site <- beta_values[rownames(beta_values) == cpg_site, , drop = FALSE]

    methylation_long <- methylation_significant_site |>
      rownames_to_column("CpG_site") |>
      pivot_longer(cols = -CpG_site,
        names_to = sample_col, values_to = "Methylation")

    melted_df <- merge(methylation_long, metadata, by = sample_col)
    cpg_metadata <- dmp_results[dmp_results[[cpg_col]] == cpg_site, ]

    plot <- ggplot(melted_df, aes(x = .data[[group_col]], y = Methylation, fill = .data[[group_col]])) +
      # geom_violin(trim = TRUE, color = "black", alpha = 0.2) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.1, height = 0, size = 2, alpha = 0.7) +
      stat_summary(fun = mean, geom = "point", shape = 23, size = 3, color = "black", fill = "red") +
      theme_minimal() +
      theme(legend.position = "none") +
      labs(title = cpg_site,
        subtitle = paste("p_adj =", format(cpg_metadata[[pvalue_col]], scientific = TRUE)),
        x = "",
        y = "Methylation beta value")

    boxplot_list[[cpg_site]] <- plot
  }

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    ggsave(file.path(output_dir, "volcano_plot.png"), volcano_plot, width = 7, height = 5, dpi = 300)

    pdf(file.path(output_dir, "dmp_boxplots.pdf"), width = 8, height = 6)
    for (plot in boxplot_list) {
      print(plot)
    }
    dev.off()
  }

  list(volcano_plot = volcano_plot, boxplots = boxplot_list)
}

run_dmr_analysis <- function(m_values,
                             metadata,
                             group_col = "stage",
                             contrast_string = "stagelate_control - stagePostConverter",
                             covariates = c("age", "sex", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"),
                             sample_col = "Sample_Name",
                             subject_col = "TREND.ID",
                             arraytype = "EPIC",
                             lambda = 1000,
                             C = 2,
                             genome = "hg19",
                             output_file = NULL) {
  install_if_missing(
    packages = character(),
    bioc_packages = c("limma", "DMRcate")
  )
  library(limma)
  library(DMRcate)

  m_values <- m_values[, colnames(m_values) %in% metadata[[sample_col]]]
  m_values <- m_values[, metadata[[sample_col]]]

  design_terms <- paste(c(paste0("0 + ", group_col), covariates), collapse = " + ")
  design <- model.matrix(as.formula(paste("~", design_terms)), data = metadata)
  contrast_matrix <- makeContrasts(contrasts = contrast_string, levels = design)

  block <- factor(metadata[[subject_col]])
  correlation <- duplicateCorrelation(as.matrix(m_values), design = design, block = block)

  cpg_annotation <- cpg.annotate(object = as.matrix(m_values), datatype = "array",
    what = "M", analysis.type = "differential",
    design = design, contrasts = TRUE,
    cont.matrix = contrast_matrix, coef = 1,
    block = block, correlation = correlation$consensus,
    arraytype = arraytype)

  dmr_results <- dmrcate(cpg_annotation, lambda = lambda, C = C)
  dmr_ranges <- extractRanges(dmr_results, genome = genome)
  dmr_table <- as.data.frame(dmr_ranges)

  if (!is.null(output_file)) {
    write.csv(dmr_table, output_file, row.names = FALSE)
  }

  list(annotation = cpg_annotation, duplicate_correlation = correlation, dmr_results = dmr_results, ranges = dmr_ranges, table = dmr_table)
}

annotate_dmrs_and_plot_heatmaps <- function(dmr_ranges,
                                            beta_values,
                                            metadata,
                                            annotation_data,
                                            sample_col = "Sample_Name",
                                            group_col = "stage",
                                            annotation_cols = c("stage", "age"),
                                            top_n = 10,
                                            output_prefix) {
  install_if_missing(
    packages = c("dplyr", "pheatmap"),
    bioc_packages = c("GenomicRanges", "IRanges")
  )
  library(dplyr)
  library(pheatmap)
  library(GenomicRanges)
  library(IRanges)

  cpg_granges <- GRanges(
    seqnames = annotation_data$chr,
    ranges = IRanges(start = annotation_data$pos, end = annotation_data$pos),
    cpg_id = rownames(annotation_data),
    gene_name = annotation_data$UCSC_RefGene_Name,
    gene_group = annotation_data$UCSC_RefGene_Group
  )

  overlaps <- findOverlaps(dmr_ranges, cpg_granges)

  overlap_df <- data.frame(
    DMR_ID = queryHits(overlaps),
    cpg_id = cpg_granges$cpg_id[subjectHits(overlaps)],
    gene_name = cpg_granges$gene_name[subjectHits(overlaps)],
    gene_group = cpg_granges$gene_group[subjectHits(overlaps)]
  )

  overlap_df$gene_name <- as.character(overlap_df$gene_name)
  overlap_df$gene_name[overlap_df$gene_name == ""] <- NA

  genes_by_dmr <- overlap_df %>%
    group_by(DMR_ID) %>%
    summarise(
      cpg.probes = paste(unique(cpg_id), collapse = ";"),
      overlapping.genes = paste(unique(na.omit(unlist(strsplit(gene_name, ";")))), collapse = ";"),
      overlapping.gene.groups = paste(unique(na.omit(unlist(strsplit(gene_group, ";")))), collapse = ";"),
      .groups = "drop"
    )

  annotated_dmrs <- as.data.frame(dmr_ranges) %>%
    mutate(DMR_ID = row_number()) %>%
    left_join(genes_by_dmr, by = "DMR_ID")

  write.csv(annotated_dmrs, paste0(output_prefix, "_annotated_dmrs.csv"), row.names = FALSE)
  write.csv(overlap_df, paste0(output_prefix, "_dmr_cpg_gene_map.csv"), row.names = FALSE)

  annotation_col <- metadata
  rownames(annotation_col) <- annotation_col[[sample_col]]
  annotation_col <- annotation_col[, annotation_cols, drop = FALSE]

  sorted_samples <- rownames(annotation_col[order(annotation_col[[group_col]]), , drop = FALSE])
  beta_values <- beta_values[, sorted_samples]

  selected_dmrs <- annotated_dmrs$DMR_ID[1:min(top_n, nrow(annotated_dmrs))]

  pdf(paste0(output_prefix, "_dmr_heatmaps.pdf"), width = 12, height = 9)

  for (dmr_id in selected_dmrs) {
    current_cpgs <- overlap_df$cpg_id[overlap_df$DMR_ID == dmr_id]
    cpgs_in_data <- intersect(current_cpgs, rownames(beta_values))

    dmr_matrix <- beta_values[cpgs_in_data, , drop = FALSE]

    gene_lookup <- overlap_df$gene_name[match(rownames(dmr_matrix), overlap_df$cpg_id)]
    gene_labels <- sapply(strsplit(gene_lookup, ";"), `[`, 1)
    rownames(dmr_matrix) <- paste0(gene_labels, "_", rownames(dmr_matrix))

    dmr_row <- annotated_dmrs[annotated_dmrs$DMR_ID == dmr_id, ]

    pheatmap(
      as.matrix(dmr_matrix),
      main = paste(
        "DMR", dmr_id,
        "\n", dmr_row$seqnames, ":", dmr_row$start, "-", dmr_row$end,
        "\nOverlapping genes:", dmr_row$overlapping.genes
      ),
      color = colorRampPalette(c("blue", "white", "red"))(100),
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      annotation_col = annotation_col,
      show_rownames = TRUE,
      show_colnames = FALSE
    )
  }

  dev.off()

  list(annotated_dmrs = annotated_dmrs, dmr_cpg_gene_map = overlap_df)
}
