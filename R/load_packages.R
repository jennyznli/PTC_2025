# load_packages.R
# Purpose: Install and load all required CRAN and Bioconductor packages

cat("Loading required packages...\n")

# ---- Setup ----
set.seed(123)
DATE <- format(Sys.Date(), "%Y%m%d")

# load_packages.R
# Purpose: Install and load all required packages

# ---- Setup ----
set.seed(123)

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", quiet = TRUE)
}

# ---- Define Required Packages ----
packages <- list(
    cran = c(
        "dplyr", "ggplot2", "readr", "tidyr", "plotly", "readxl", "here", "stringr",
        "Rtsne", "pheatmap", "RColorBrewer", "ggpubr", "ggrepel", "clusterProfiler",
        "msigdbr", "fgsea", "enrichplot", "pvclust", "tidyverse", "randomForest",
        "sail", "RSNNS", "e1071", "caret", "ISLR", "pROC", "shapviz", "future",
        "future.apply", "logger", "treeshap", "viridis", "pals"
    ),
    bioc = c(
        "BiocParallel", "sesame", "DESeq2", "limma", "AnnotationHub", "CytoMethIC",
        "SummarizedExperiment", "knowYourCG"
    )
)

# ---- Install and Load Function ----
load_packages <- function(pkg_list) {
    # Install missing CRAN packages
    missing_cran <- pkg_list$cran[!sapply(pkg_list$cran, requireNamespace, quietly = TRUE)]
    if (length(missing_cran) > 0) {
        install.packages(missing_cran, quiet = TRUE)
    }

    # Install missing Bioconductor packages
    missing_bioc <- pkg_list$bioc[!sapply(pkg_list$bioc, requireNamespace, quietly = TRUE)]
    if (length(missing_bioc) > 0) {
        BiocManager::install(missing_bioc, ask = FALSE, update = FALSE, quiet = TRUE)
    }

    # Load all packages
    all_pkgs <- c(pkg_list$cran, pkg_list$bioc)
    invisible(sapply(all_pkgs, function(x) {
        suppressPackageStartupMessages(library(x, character.only = TRUE))
    }))

    message("All packages loaded successfully!")
}

# ---- Execute ----
load_packages(packages)

# ---- Confirm 'here' is working ----
here::here()
