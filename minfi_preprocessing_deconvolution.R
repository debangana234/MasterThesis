source("package_helpers.R")

preprocess_minfi_data <- function(data_dir,
                                  output_dir,
                                  sample_sheet = "completeSample_sheet.csv",
                                  detection_p_cutoff = 0.01) {
  install_if_missing(packages = c("RColorBrewer"), bioc_packages = c("minfi")
  )
  library(minfi)
  library(RColorBrewer)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  targets <- read.metharray.sheet(data_dir, pattern = sample_sheet)
  rgset <- read.metharray.exp(targets = targets)

  mset_raw <- preprocessRaw(rgset)
  rset_raw <- ratioConvert(mset_raw, what = "both", keepCN = TRUE)
  gset_raw <- mapToGenome(rset_raw)

  raw_beta_values <- as.data.frame(getBeta(gset_raw))
  raw_m_values <- as.data.frame(getM(gset_raw))
  copy_number_values <- as.data.frame(getCN(gset_raw))

  det_p <- detectionP(rgset)
  mean_detection_p <- colMeans(det_p)

  pdf(file.path(output_dir, "minfi_qc_detection_p.pdf"), width = 8, height = 6)
  plotQC(getQC(mset_raw))
  barplot(mean_detection_p, las = 2, cex.names = 0.8, ylab = "Mean detection p-values")
  abline(h = 0.05, col = "red")
  dev.off()

  mset_norm <- preprocessQuantile(
    rgset,
    fixOutliers = TRUE,
    removeBadSamples = TRUE,
    badSampleCutoff = 10.5,
    quantileNormalize = TRUE,
    stratified = TRUE,
    mergeManifest = FALSE,
    sex = NULL
  )

  beta_values <- as.data.frame(getBeta(mset_norm))
  m_values <- as.data.frame(getM(mset_norm))

  colnames(beta_values) <- targets$Sample_Name
  colnames(m_values) <- targets$Sample_Name
  colnames(raw_beta_values) <- targets$Sample_Name
  colnames(raw_m_values) <- targets$Sample_Name
  colnames(copy_number_values) <- targets$Sample_Name

  keep_probes <- rowSums(det_p < detection_p_cutoff) == ncol(mset_norm)
  mset_filtered <- mset_norm[keep_probes, ]

  mset_no_snps <- dropLociWithSnps(mset_filtered, snps = c("SBE", "CpG"), maf = 0)
  annotation_data <- getAnnotation(mset_no_snps)

  sex_chr_probes <- annotation_data[annotation_data$chr %in% c("chrX", "chrY"), ]
  sex_chr_probe_ids <- rownames(sex_chr_probes)

  beta_no_snps <- getBeta(mset_no_snps)
  m_values_no_snps <- getM(mset_no_snps)
  probes_to_keep <- setdiff(rownames(beta_no_snps), sex_chr_probe_ids)

  beta_no_sex <- beta_no_snps[probes_to_keep, ]
  m_values_no_sex <- m_values_no_snps[probes_to_keep, ]

  cpg_probes <- annotation_data[annotation_data$Type %in% c("I", "II"), ]
  cpg_probe_names <- cpg_probes$Name

  beta_filtered <- as.data.frame(getBeta(mset_filtered))
  m_values_filtered <- as.data.frame(getM(mset_filtered))
  beta_clean <- as.data.frame(beta_no_sex[rownames(beta_no_sex) %in% cpg_probe_names, ])
  m_values_clean <- as.data.frame(m_values_no_sex[rownames(m_values_no_sex) %in% cpg_probe_names, ])

  colnames(beta_filtered) <- targets$Sample_Name
  colnames(m_values_filtered) <- targets$Sample_Name
  colnames(beta_clean) <- targets$Sample_Name
  colnames(m_values_clean) <- targets$Sample_Name

  preprocessing_summary <- data.frame(
    step = c(
      "raw_rgset",
      "quantile_normalised",
      "detection_p_filtered",
      "snp_filtered",
      "sex_chr_and_cpg_filtered"
    ),
    probes = c(nrow(rgset), nrow(mset_norm),
      nrow(mset_filtered), nrow(mset_no_snps),
      nrow(beta_clean)),
    samples = c(
      ncol(rgset),
      ncol(mset_norm),
      ncol(mset_filtered),
      ncol(mset_no_snps),
      ncol(beta_clean)
    )
  )

  write.csv(beta_values, file.path(output_dir, "beta_quantile_normalised.csv"))
  write.csv(beta_filtered, file.path(output_dir, "beta_detection_p_filtered.csv"))
  write.csv(beta_clean, file.path(output_dir, "beta_clean.csv"))
  write.csv(m_values, file.path(output_dir, "m_values_quantile_normalised.csv"))
  write.csv(m_values_clean, file.path(output_dir, "m_values_clean.csv"))
  write.csv(preprocessing_summary, file.path(output_dir, "preprocessing_summary.csv"), row.names = FALSE)
  write.csv(data.frame(Sample_Name = names(mean_detection_p), mean_detection_p = mean_detection_p),
            file.path(output_dir, "mean_detection_p.csv"),
            row.names = FALSE)

  saveRDS(rgset, file.path(output_dir, "rgset.rds"))
  saveRDS(mset_raw, file.path(output_dir, "mset_raw.rds"))
  saveRDS(mset_norm, file.path(output_dir, "mset_norm.rds"))
  saveRDS(mset_filtered, file.path(output_dir, "mset_filtered_detection_p.rds"))
  saveRDS(mset_no_snps, file.path(output_dir, "mset_no_snps.rds"))
  saveRDS(beta_clean, file.path(output_dir, "beta_clean.rds"))
  saveRDS(m_values_clean, file.path(output_dir, "m_values_clean.rds"))

  list(
    targets = targets,
    summary = preprocessing_summary,
    mean_detection_p = mean_detection_p,
    beta_clean = beta_clean,
    m_values_clean = m_values_clean
  )
}

deconvolute_epic_reference <- function(rgset, sample_names,
                                       cell_types = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"),
                                       output_file) {
  install_if_missing(
    packages = character(),
    bioc_packages = c("minfi", "FlowSorted.Blood.EPIC")
  )
  library(minfi)
  library(FlowSorted.Blood.EPIC)

  cell_counts <- estimateCellCounts2(
    rgset,
    compositeCellType = "Blood",
    processMethod = "preprocessNoob",
    probeSelect = "IDOL",
    cellTypes = cell_types,
    referencePlatform = "IlluminaHumanMethylationEPIC",
    referenceset = NULL,
    CustomCpGs = NULL,
    returnAll = FALSE,
    meanPlot = FALSE,
    verbose = TRUE
  )

  cell_counts <- as.data.frame(cell_counts$prop)
  rownames(cell_counts) <- sample_names

  write.csv(cell_counts, output_file)

  cell_counts
}
