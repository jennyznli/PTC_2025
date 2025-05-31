# ============================================================
#   This script performs preprocessing for
#   pediatric + adult betas.
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
# LOADING DATA
# ========================
betas_ped <- readRDS(file.path(DATA_DIR, "ped_betas_processed.rds"))
cat("Original uncollapsed pediatric EPICv2 beta dimensions:", dim(betas_ped), "\n")
betas_ped <- betasCollapseToPfx(betas_ped, BPPARAM = BiocParallel::MulticoreParam(N_WORKERS))
cat("Collapsed EPICv2 beta dimensions:", dim(betas_ped), "\n")

betas_adult <- readRDS(file.path(DATA_DIR, "adult_betas_processed.rds"))
cat("Original HM450 beta dimensions:", dim(betas_adult), "\n")

ss_ped <- read_excel(file.path(SS_DIR, "pediatric_master.xlsx"))
ss_adult <- read_excel(file.path(SS_DIR, "adult_master.xlsx"))

# ========================
# COMBINING BETAS
# ========================
ped_probes <- rownames(betas_ped)
adult_probes <- rownames(betas_adult)
cat("Pediatric probe count:", length(ped_probes), "\n")
cat("Adult probe count:", length(adult_probes), "\n")

joint_probes <- intersect(ped_probes, adult_probes)
cat("Overlapping probes:", length(joint_probes), "\n")

joint_betas <- cbind(betas_ped[joint_probes,], betas_adult[joint_probes,])
cat("Joint beta matrix dimensions:", dim(joint_betas), "\n")

saveRDS(joint_betas, file.path(DATA_DIR, "joint_betas_processed.rds"))

# ========================
# IMPUTATION
# ========================
# clean and impute
joint_betas_clean <- cleanMatrixForClusterW(joint_betas)
joint_betas_imputed <- impute(joint_betas_clean, "HM450")

# save primary + tumor
cat("Joint beta matrix dimensions:", dim(joint_betas_imputed), "\n")
saveRDS(joint_betas_imputed, file.path(DATA_DIR, "joint_betas_imputed.rds"))
