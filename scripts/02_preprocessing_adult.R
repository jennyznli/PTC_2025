# ============================================================
#   This script performs preprocessing for
#   adult betas.
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
ss = read_excel(file.path(SS_DIR, "adult_master.xlsx"))

thca = load(file.path(DATA_DIR, "THCA.rda"))
betas <- get("betas")

# ========================
# PREPROCESSING
# ========================
colnames(betas) <- substr(colnames(betas), 1, 15)
betas = betas[, colnames(betas) %in% ss$Sample_ID]

# get masks
masks <- getMask(platform = "HM450", mask_names = "recommended")
cat("Found", length(masks), "HM450 probes to be masked.\n")

# get sex chromosome probes
mft <- sesameAnno_readManifestTSV("HM450.hg38.manifest")
sex <- mft %>%
    filter(CpG_chrm %in% c("chrX", "chrY"))
sex_probes <- sex$Probe_ID
cat("Found", length(sex_probes), "HM450 sex probes.\n")

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

betas <- betas[, match(ss$Sample_ID, colnames(betas))]
saveRDS(betas, file.path(DATA_DIR, "adult_betas_processed.rds"))
cat("After preprocessing beta dimensions:", dim(betas), "\n")

# ========================
# IMPUTATION
# ========================
# clean and impute
betas_clean <- cleanMatrixForClusterW(betas)
betas_imputed <- impute(betas_clean, "HM450")
colnames(betas_imputed) <- ss$Sample_ID

# save primary + tumor
cat("Final beta dimension:", dim(betas_imputed), "\n")
saveRDS(betas_imputed, file.path(DATA_DIR, "adult_betas_imputed.rds"))
