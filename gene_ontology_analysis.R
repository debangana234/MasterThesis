source("Parkinson/clean_scripts/00_package_helpers.R")

run_go_profiler_top500 <- function(top500_cpgs,
                                   gene_col = "UCSC_RefGene_Name",
                                   organism = "hsapiens",
                                   domain_scope = "annotated",
                                   correction_method = "g_SCS",
                                   output_file) {
  install_if_missing(packages = c("gprofiler2", "dplyr"))
  library(gprofiler2)
  library(dplyr)

  genes <- top500_cpgs[[gene_col]]
  genes <- genes[genes != ""]
  genes <- na.omit(genes)
  genes <- unlist(strsplit(genes, ";"))
  genes <- trimws(genes)
  genes <- unique(genes)
  genes <- genes[genes != ""]

  gost_results <- gost(
    query = genes,
    organism = organism,
    ordered_query = FALSE,
    multi_query = FALSE,
    significant = TRUE,
    domain_scope = domain_scope,
    correction_method = correction_method
  )

  results_df <- gost_results$result
  results_df <- results_df %>%
    mutate(across(where(is.list), ~ sapply(.x, toString)))

  write.csv(results_df, output_file, row.names = FALSE)

  list(genes = genes, gost_results = gost_results, results = results_df)
}

map_go_genes_to_cpgs <- function(top500_cpgs,
                                 go_results,
                                 cpg_col = "cpg_id",
                                 gene_col = "UCSC_RefGene_Name",
                                 output_file) {
  install_if_missing(packages = c("dplyr"))
  library(dplyr)

  go_genes <- go_results$intersection
  go_genes <- as.character(go_genes)
  go_genes <- unlist(strsplit(go_genes, ","))
  go_genes <- trimws(go_genes)
  go_genes <- unique(go_genes)
  go_genes <- go_genes[go_genes != ""]

  top500_cpgs[[gene_col]] <- as.character(top500_cpgs[[gene_col]])
  cpg_gene_list <- strsplit(top500_cpgs[[gene_col]], ";")

  keep_cpgs <- sapply(cpg_gene_list, function(genes_for_cpg) {
    genes_for_cpg <- trimws(genes_for_cpg)
    any(genes_for_cpg %in% go_genes)
  })

  mapped_cpgs <- top500_cpgs[keep_cpgs, ]
  mapped_gene_list <- strsplit(mapped_cpgs[[gene_col]], ";")

  mapped_cpgs$mapped_go_genes <- sapply(mapped_gene_list, function(genes_for_cpg) {
    genes_for_cpg <- trimws(genes_for_cpg)
    paste(intersect(genes_for_cpg, go_genes), collapse = ";")
  })

  write.csv(mapped_cpgs, output_file, row.names = FALSE)

  list(go_genes = go_genes, mapped_cpgs = mapped_cpgs)
}

plot_go_bubbleplot <- function(go_results, output_prefix, top_n = 20) {
  install_if_missing(packages = c("dplyr", "ggplot2", "stringr"))
  library(dplyr)
  library(ggplot2)
  library(stringr)

  plot_data <- go_results %>%
    arrange(p_value) %>%
    head(top_n) %>%
    mutate(log10_adjusted_p = -log10(p_value),
      term_name_clean = str_wrap(term_name, width = 34)) %>%
    arrange(log10_adjusted_p)

  plot_data$term_name_clean <- factor(plot_data$term_name_clean, levels = plot_data$term_name_clean)

  go_plot <- ggplot(plot_data, aes(x = log10_adjusted_p, y = term_name_clean)) +
    geom_point(aes(size = intersection_size, fill = source),
      shape = 21, color = "black",
      alpha = 0.9, stroke = 0.4) +
    scale_size_continuous(range = c(3, 11)) +
    labs(x = "-log10 adjusted p-value",
      y = "", size = "Genes in term", fill = "Source") +
    theme_minimal() +
    theme(axis.text.y = element_text(color = "black"),
      panel.grid.minor = element_blank(),
      legend.position = "right")

  write.csv(plot_data, paste0(output_prefix, "_bubbleplot_data.csv"), row.names = FALSE)
  ggsave(paste0(output_prefix, ".png"), go_plot, width = 10, height = 7, dpi = 300)
  ggsave(paste0(output_prefix, ".pdf"), go_plot, width = 10, height = 7)

  go_plot
}
