# ============================================================
# Install and load all required CRAN and Bioconductor packages
# ============================================================

cat("Loading required packages...\n")

set.seed(123)
DATE <- format(Sys.Date(), "%Y%m%d")

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", quiet = TRUE)
}

packages <- list(
    cran = c(
        "dplyr", "ggplot2", "readr", "tidyr", "plotly", "readxl", "here", "stringr",
        "Rtsne", "pheatmap", "RColorBrewer", "ggpubr", "ggrepel", "clusterProfiler",
        "msigdbr", "fgsea", "enrichplot", "pvclust", "tidyverse", "randomForest",
        "RSNNS", "e1071", "caret", "ISLR", "pROC", "shapviz", "future",
        "future.apply", "logger", "treeshap", "viridis", "pals", "openxlsx"
    ),
    bioc = c(
        "BiocParallel", "sesame", "DESeq2", "limma", "AnnotationHub", "CytoMethIC",
        "SummarizedExperiment", "knowYourCG"
    )
)

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

load_packages(packages)

