#!/usr/bin/env Rscript

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cran.r-project.org")
}

# Install Bioconductor 3.18 (compatible with R 4.3)
BiocManager::install(version = "3.18")

# Install GenomicRanges (will get the compatible version for R 4.3)
BiocManager::install("GenomicRanges")
