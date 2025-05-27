# ============================================================
#   This script performs preprocessing for
#   pediatric betas.
# ============================================================

# ========================
# CONFIGURATION
# ========================
N_WORKERS <- 20

BASE_DIR <- "/home/lijz/thyroid"
R_DIR <- file.path(BASE_DIR, "R")
SS_DIR <- file.path(BASE_DIR, "ss")
DATA_DIR <- file.path(BASE_DIR, "data")
IDAT_DIR <- file.path(BASE_DIR, "IDAT")

source(file.path(R_DIR, "functions.R"))
source(file.path(R_DIR, "color_keys.R"))
source(file.path(R_DIR, "load_packages.R"))

# ========================
# LOADING DATA
# ========================
ss <- read_excel(file.path(SS_DIR, "pediatric_master.xlsx"))
ss_primary <- ss %>% filter(Lymph_Node == "F")

# ========================
# PREPROCESSING
# ========================
betas <- openSesame(IDAT_DIR, prep = "QCDPB", func = getBetas, BPPARAM = BiocParallel::MulticoreParam(N_WORKERS))
cat("Raw beta dimensions:", dim(betas), "\n")
saveRDS(betas, file.path(DATA_DIR, "ped_betas_unprocessed.rds"))

# applying additional masks
masks <- getMask(platform = "EPICv2", mask_names = "recommended")
cat("Found", length(masks), "EPICv2 probes to be masked.\n")

# get sex chromosome probes
mft <- sesameAnno_readManifestTSV("EPICv2.hg38.manifest")
sex <- mft %>%
    filter(CpG_chrm %in% c("chrX", "chrY"))
sex_probes <- sex$Probe_ID
cat("Found", length(sex_probes), "EPICv2 sex probes.\n")

# remove masked and sex chromosome probes
remove_probes <- unique(c(sex_probes, masks))
betas <- betas[!rownames(betas) %in% remove_probes, ]
cat("After removing sex and masked probes beta dimensions:", dim(betas), "\n")

# probe types
probes <- rownames(betas)
probe_types <- table(substr(probes, 1, 2))
cat("Probe types before filtering:\n")
print(probe_types)

# filter to only cg probes
cg_probes <- probes[substr(probes, 1, 2) == "cg"]
betas <- betas[cg_probes, ]
cat("After keeping only CG probes beta dimensions:", dim(betas), "\n")
cat("Kept", length(cg_probes), "CG probes out of", length(probes), "total probes\n")

betas <- betas[, match(ss$IDAT, colnames(betas))]
saveRDS(betas, file.path(DATA_DIR, "ped_betas_processed.rds"))
cat("After preprocessing beta dimensions:", dim(betas), "\n")

# ========================
# IMPUTATION
# ========================
betas_clean <- cleanMatrixForClusterW(betas)
betas_imputed <- impute(betas_clean, "EPICv2")

# collapse to pfx
betas_collapsed <- betasCollapseToPfx(betas_imputed, BPPARAM = BiocParallel::MulticoreParam(N_WORKERS))
cat("After collapsing to prefixes:", dim(betas_collapsed), "\n")

# save primary + tumor
cat("Final primary & LN beta dimension:", dim(betas_collapsed), "\n")
saveRDS(betas_collapsed, file.path(DATA_DIR, "ped_betas_imputed.rds"))

# save primary only
betas_primary <- betas_collapsed[, match(ss_primary$IDAT, colnames(betas_collapsed))]
cat("Final primary only beta dimension:", dim(betas_primary), "\n")

saveRDS(betas_primary, file.path(DATA_DIR, "ped_betas_primary_imputed.rds"))

