# ============================================================
#   FIG 2.
#   This script performs differential methylation analysis
#   for 88 primary tumor samples:
#   Invasiveness & Methylation Clusters
# ============================================================
# ========================
# CONFIGURATION
# ========================
N_WORKERS <- 20

BASE_DIR <- "/home/lijz/thyroid"
R_DIR <- file.path(BASE_DIR, "R")
SS_DIR <- file.path(BASE_DIR, "ss")
DATA_DIR <- file.path(BASE_DIR, "data")

source(file.path(R_DIR, "functions.R"))
source(file.path(R_DIR, "color_keys.R"))
source(file.path(R_DIR, "load_packages.R"))

# ========================
# INVASIVENESS
# ========================
cat("=== CLINICAL INVASIVENESS ANALYSIS ===\n")

# Load data
ss <- read_excel(file.path(SS_DIR, "pediatric_master.xlsx"))
ss_primary <- ss %>% filter(Lymph_Node == "F")
betas <- readRDS(file.path(DATA_DIR, "ped_betas_primary_imputed.rds"))

ss_primary$Age <- as.numeric(ss_primary$Age)
ss_primary$Sex <- as.factor(ss_primary$Sex)
ss_primary$Clinical_Invasiveness <- as.factor(ss_primary$Clinical_Invasiveness)

se <- SummarizedExperiment(betas, colData = ss_primary)

se_ok <- (checkLevels(assay(se), colData(se)$Sex) &
              checkLevels(assay(se), colData(se)$Clinical_Invasiveness))

cat("Probes passing QC:", sum(se_ok), "out of", length(se_ok), "\n")
se <- se[se_ok, ]
cat("Final SE dimensions:", dim(se), "\n")

colData(se)$Sex <- relevel(factor(colData(se)$Sex), "Female")
colData(se)$Clinical_Invasiveness <- relevel(factor(colData(se)$Clinical_Invasiveness), "Low")

cat("Reference levels - Sex:", levels(colData(se)$Sex)[1],
    "Clinical_Invasiveness:", levels(colData(se)$Clinical_Invasiveness)[1], "\n")

smry <- DML(se, ~Clinical_Invasiveness + Sex + Age,
            BPPARAM = BiocParallel::MulticoreParam(workers = N_WORKERS))
res_inv <- summaryExtractTest(smry)

cat("Results dimensions:", dim(res_inv), "\n")
cat("First few results:\n")
print(head(res_inv, 3))

saveRDS(res_inv, file = file.path(DATA_DIR, "ped88_dm_invasiveness.rds"))
cat("Clinical_Invasiveness differential methylation completed:", dim(res_inv), "\n")

# ========================
# CLUSTERS: LI REFERENCE
# ========================
cat("=== LI REFERENCE CLUSTER ANALYSIS ===\n")

ss <- read_excel(file.path(SS_DIR, "pediatric_master.xlsx"))
ss_primary <- ss %>% filter(Lymph_Node == "F")
betas <- readRDS(file.path(DATA_DIR, "ped_betas_primary_imputed.rds"))

ss_primary$Age <- as.numeric(ss_primary$Age)
ss_primary$Sex <- as.factor(ss_primary$Sex)
ss_primary$Methylation_Clusters <- as.factor(ss_primary$Methylation_Clusters)

cat("Sample groups:", table(ss_primary$Methylation_Clusters), "\n")

se <- SummarizedExperiment(betas, colData = ss_primary)

se_ok <- (checkLevels(assay(se), colData(se)$Sex) &
              checkLevels(assay(se), colData(se)$Methylation_Clusters))

cat("Probes passing QC:", sum(se_ok), "out of", length(se_ok), "\n")
se <- se[se_ok, ]
cat("Final SE dimensions:", dim(se), "\n")

colData(se)$Sex <- relevel(factor(colData(se)$Sex), "Female")
colData(se)$Methylation_Clusters <- relevel(factor(colData(se)$Methylation_Clusters), "LI")

cat("LI reference levels - Sex:", levels(colData(se)$Sex)[1],
    "Methylation_Clusters:", levels(colData(se)$Methylation_Clusters)[1], "\n")

smry <- DML(se, ~Methylation_Clusters + Sex + Age,
            BPPARAM = BiocParallel::MulticoreParam(workers = N_WORKERS))
res_li <- summaryExtractTest(smry)

cat("Results dimensions:", dim(res_li), "\n")
cat("First few results:\n")
print(head(res_li, 3))

saveRDS(res_li, file = file.path(DATA_DIR, "ped88_dm_cluster_li.rds"))
cat("Cluster (LI reference) differential methylation completed:", dim(res_li), "\n")

# ========================
# CLUSTERS: HIL REFERENCE
# ========================
cat("=== HIL REFERENCE CLUSTER ANALYSIS ===\n")

ss <- read_excel(file.path(SS_DIR, "pediatric_master.xlsx"))
ss_primary <- ss %>% filter(Lymph_Node == "F")
betas <- readRDS(file.path(DATA_DIR, "ped_betas_primary_imputed.rds"))

ss_primary$Age <- as.numeric(ss_primary$Age)
ss_primary$Sex <- as.factor(ss_primary$Sex)
ss_primary$Methylation_Clusters <- as.factor(ss_primary$Methylation_Clusters)

se <- SummarizedExperiment(betas, colData = ss_primary)

se_ok <- (checkLevels(assay(se), colData(se)$Sex) &
              checkLevels(assay(se), colData(se)$Methylation_Clusters))

cat("Probes passing QC:", sum(se_ok), "out of", length(se_ok), "\n")
se <- se[se_ok, ]
cat("Final SE dimensions:", dim(se), "\n")

colData(se)$Sex <- relevel(factor(colData(se)$Sex), "Female")
colData(se)$Methylation_Clusters <- relevel(factor(colData(se)$Methylation_Clusters), "HIL")

cat("HIL reference levels - Sex:", levels(colData(se)$Sex)[1],
    "Methylation_Clusters:", levels(colData(se)$Methylation_Clusters)[1], "\n")

cat("Starting HIL reference DML analysis...\n")
smry <- DML(se, ~Methylation_Clusters + Sex + Age,
            BPPARAM = BiocParallel::MulticoreParam(workers = N_WORKERS))
res_hil <- summaryExtractTest(smry)

cat("Results dimensions:", dim(res_hil), "\n")
cat("First few results:\n")
print(head(res_hil, 3))

saveRDS(res_hil, file = file.path(DATA_DIR, "ped88_dm_cluster_hil.rds"))
cat("Cluster (HIL reference) differential methylation completed:", dim(res_hil), "\n")

# ========================
# MULTIPLE TESTING CORRECTION
# ========================
cat("=== MULTIPLE TESTING CORRECTION ===\n")

cat("Processing invasiveness results for BH correction...\n")
if ("Pval_Clinical_InvasivenessHigh" %in% colnames(res_inv)) {
    res_inv$BH_Clinical_InvasivenessHigh <- p.adjust(res_inv$Pval_Clinical_InvasivenessHigh, method = "BH")
    cat("Adjusted", nrow(res_inv), "p-values for Clinical_Invasiveness\n")
    saveRDS(res_inv, file.path(DATA_DIR, "ped88_dm_invasiveness_bh.rds"))
} else {
    cat("Warning: Pval_Clinical_InvasivenessHigh column not found in Clinical_Invasiveness results\n")
}

# Cluster groups BH correction (across all comparisons)
cat("Processing cluster group results for joint BH correction...\n")
# Combine all p-values for joint correction
all_pvals <- data.frame(
    cpg = c(rownames(res_li), rownames(res_li), rownames(res_hil)),
    comparison = c(
        rep("LIvHI", nrow(res_li)),
        rep("LIvHIL", nrow(res_li)),
        rep("HILvHI", nrow(res_hil))
    ),
    pval = c(
        res_li$Pval_Methylation_ClustersHI,    # HI vs LI
        res_li$Pval_Methylation_ClustersHIL,   # HIL vs LI
        res_hil$Pval_Methylation_ClustersHI    # HI vs HIL
    )
)

cat("Total p-values for joint correction:", nrow(all_pvals), "\n")

# Apply joint BH correction
adjusted_all_pvals <- p.adjust(all_pvals$pval, method = "BH")

# Distribute corrected p-values back to results
n <- nrow(res_li)
res_li$BH_Methylation_ClustersHI <- adjusted_all_pvals[1:n]
res_li$BH_Methylation_ClustersHIL <- adjusted_all_pvals[(n+1):(2*n)]
res_hil$BH_Methylation_ClustersHI <- adjusted_all_pvals[(2*n+1):(3*n)]

saveRDS(res_li, file.path(DATA_DIR, "ped88_dm_cluster_li_bh.rds"))
saveRDS(res_hil, file.path(DATA_DIR, "ped88_dm_cluster_hil_bh.rds"))
