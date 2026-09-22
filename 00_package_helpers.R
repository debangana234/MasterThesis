install_if_missing <- function(packages, bioc_packages = character()) {
  missing_cran <- packages[!packages %in% rownames(installed.packages())]

  if (length(missing_cran) > 0) {
    install.packages(missing_cran)
  }

  missing_bioc <- bioc_packages[!bioc_packages %in% rownames(installed.packages())]

  if (length(missing_bioc) > 0) {
    if (!"BiocManager" %in% rownames(installed.packages())) {
      install.packages("BiocManager")
    }
    library(BiocManager)
    install(missing_bioc)
  }
}
