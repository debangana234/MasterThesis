
This repository contains the R scripts used for DNA methylation analysis in my Master's thesis project on Parkinson's disease progression.

The workflow includes methylation preprocessing, immune cell type deconvolution, immune proportion analysis, epigenetic age analysis, differential methylation analysis, gene ontology analysis, and predictive modelling.

Immune cell type deconvolution was performed using the EPIC reference dataset. In addition to the TREND dataset used for the main analysis, two public Parkinson's disease methylation datasets were also used for comparison and validation: GSE111629 and GSE145361.

## Repository Structure

| Script | Description |
|---|---|
| `package_helpers.R` | Helper function for installing and loading required R packages. |
| `minfi_preprocessing_deconvolution.R` | Minfi-based preprocessing of Illumina methylation array data and reference-based blood cell type deconvolution using the EPIC reference dataset. |
| `immune_proportions_analysis.R` | Analysis of estimated immune cell type proportions, including boxplots, mixed models, and diagnosis-aligned immune proportion changes. |
| `epigenetic_age_analysis.R` | Application of epigenetic clocks, epigenetic age acceleration analysis, age acceleration boxplots, and PCA biplot. |
| `differential_methylation_analysis.R` | Differential methylation analysis using limma, including adjustment for age, sex, and estimated immune cell proportions. |
| `dmp_dmr_plots.R` | Volcano plots, CpG-level methylation boxplots, DMR analysis, DMR annotation, and DMR heatmaps. |
| `gene_ontology_analysis.R` | Gene ontology enrichment analysis using g:Profiler |
| `predictive_modelling.R` | Predictive modelling using logistic classification and SVM models, including ROC-AUC and confusion matrix. |

## Data Availability

The raw methylation data, processed methylation matrices, metadata files, and generated output files are not included in this repository because they may contain sensitive or large data.

For executing these scripts please provide:

- Illumina methylation array IDAT files
- sample metadata
- processed beta or M-value matrices where required
- immune cell type proportion estimates where required

The public datasets referred to in the analysis are:

- GSE111629
- GSE145361

Reference-based immune cell deconvolution was carried out using the EPIC blood reference dataset through the `estimateCellCounts2` method.

## Main Analysis Steps

1. Preprocess methylation array data using `minfi`.
2. Estimate blood immune cell type proportions using the EPIC reference dataset.
3. Analyse immune cell proportions across disease stages and relative to diagnosis time.
4. Estimate epigenetic age using Horvath and BLUP clocks.
5. Calculate epigenetic age acceleration.
6. Identify differentially methylated positions using `limma`.
7. Identify and annotate differentially methylated regions using `DMRcate`.
8. Perform gene ontology enrichment analysis with g:Profiler.
9. Build predictive models using selected CpGs and immune cell proportions.
